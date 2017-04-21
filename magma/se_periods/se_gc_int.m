/******************************************************************************

 Copyright (C) 2017 Christian Neurohr

 ******************************************************************************/

// Import functions
import "superelliptic.m": RS_SEIntegrate, RS_SEIntegrate3, RS_NRootAC, RS_SEMST, RS_SEMST2, RS_MakeCCVectors, RS_SETau3, RS_SEInfTau;
import "se_anal_cont.m": AC_mthRoot;
import "se_help_funcs.m": MakeCCVector, PolynomialShiftVector;

// Low precision complex field
C20<i> := ComplexField(20);
Pi20 := Real(Pi(C20));

////////////////////////////////////////////////////////////
// ***** Parameters for Gauss-Chebychev integration ***** //
////////////////////////////////////////////////////////////

function UpperBound_M(M_0,r,V_r,Len)
	M_r := 1;
	K := 0;
	for k in [1..Len] do
		if Abs(r-V_r[k]) lt 10^-10 then
			K +:= 1;
		else 
			M_r *:= V_r[k] - r;
		end if;
	end for;
	M_r := M_0 - (1/2) * Log(M_r);
	return K,M_r;
end function;

function GC_Params(CCV,Len,n,Prec)
// Computes n and e
	// Compute r_k and r_0
	V_r := []; r_0 := 4.;
	ChangeUniverse(~CCV,C20);
	for k in [1..Len] do
		r_k := (1/2) * ( Abs(CCV[k]+1) + Abs(CCV[k]-1) );
		assert r_k gt 1;
		r_0 := Min(r_0,r_k);
		Append(~V_r,r_k);
	end for;

	// M_0 = D + Log(2Pi) + n*Log(r_0)
	M_0 := Log(2*Pi20) + n*Log(r_0) + Log(10)*Prec;

   	// Choose r
	r := (1/2) * (r_0 + 1);

	K,M_r := UpperBound_M(M_0,r,V_r,Len);
	while K ne 0 do
		A := 1 + 2*M_r/K;
		r := Cosh(r_0 - r_0/(Log(A/Argcosh(r_0))+A));
		K,M_r := UpperBound_M(M_0,r,V_r,Len);
	end while;

	// Number of integration points
	//r := Argcosh(r); // r' = Log(r+Sqrt(r^2-1))
	//M_r := (M_r/(2*r));
	//N := Ceiling(M_r);
	N := Ceiling(M_r/(2*Argcosh(r)));

	// Error
	//e :=  Exp(M_r)/(Exp(2*N*r)-1);

	return N;
end function;

function GC_Params_Tree( STree,Prec);
	return GC_Params(STree`Data[STree`WorstEdge],STree`Len-1,STree`Len+1,Prec);
end function;


/////////////////////////////////////////////
// ***** Gauss-Chebychev integration ***** //
/////////////////////////////////////////////


// CurveData = [n,g,lc_m,C]
function GC_Intergrals_Factor( VectorIntegral, ba2, cab, up, CD )
// Multiply vector integral by constants
	
	// Polynomial shift by (b+a)/(b-a)
	PolynomialShiftVector( ~VectorIntegral, cab, CD[2], 0 );
	
	// c1 = lc_m * 2 / ((b+a)/2)^(n/2)
	n2 := Floor(CD[1]/2);
	c1 := ba2^n2;
	if CD[1] mod 2 eq 1 then
		t := Sqrt(ba2);
		c1 *:= t;
	end if;
	c1 := CD[3]/c1;
	// * -I
	if up mod 2 eq 0 then
		c1 *:= -CD[4].1;
	end if;

	VectorIntegral := [ c1 * v : v in VectorIntegral ];
	// c2 = ((b+a)/2)^(i+1)
	c2 := 1;
	for k in [1..CD[2]] do
		c2 *:= ba2;
		VectorIntegral[k] *:= c2;
	end for;
	ElemInts := [ (1/2) * v : v in VectorIntegral ];
	return VectorIntegral, ElemInts;
end function;


function GC_Integrals( CCV, up, g, n, C, N)
// Gauss-Chebychev integration
	VecInt := [ Zero(C) : j in [1..g] ]; i := C.1; PiR := Real(Pi(C)); const := 1/(2*N)*PiR;
	Zetas := [1,i,-1,-i];
	for k in [1..N] do
		// Compute x
		x := Cos((2*k-1)*const);

		// Compute 1/y(x)
		y := 1/AC_mthRoot(x,CCV,Zetas,up,2,n);

		// Compute all differentials
		yxi := y;
		VecInt[1] +:= yxi;
		for j in [2..g] do
			yxi *:= x;
			VecInt[j] +:= yxi;
		end for;
	end for;

	// Multiply by weight Pi/n
	w := PiR/N;
	for j in [1..g] do
		VecInt[j] *:= w;
	end for;
	return VecInt;
end function;


// CurveData = [n,g,lc_m,C]

function GC_Integrals_Edge( EdgeData, CurveData, N)
// Integrate an edge (of the spanning tree) with Gauss-Chebychev integration }

	// Compute integral on [-1,1]
	VectorIntegral := GC_Integrals( EdgeData, Floor(EdgeData[CurveData[1]+1]), CurveData[2], CurveData[1]-2, CurveData[4], N);

	// Multiply with constants
	PeriodsEdge, ElemIntegralEdge := GC_Intergrals_Factor( VectorIntegral, EdgeData[CurveData[1]-1], EdgeData[CurveData[1]], Floor(EdgeData[CurveData[1]+1]), CurveData);

	return PeriodsEdge, ElemIntegralEdge; 
end function;


function GC_Integrals_Tree(SEC)
// Compute integrals for tree
	// Compute integration parameters for spanning tree
	N := GC_Params_Tree(SEC`SpanningTree,SEC`Prec);
	
	CurveData := <SEC`Degree[2],SEC`Genus,1,SEC`ComplexField>;
	Periods := []; ElementaryIntegrals := [];
	for k in [1..SEC`SpanningTree`Len] do
		Period, EI := GC_Integrals_Edge(SEC`SpanningTree`Data[k],CurveData,N);
		Append(~Periods,[ [ P ] : P in Period ]);
		Append(~ElementaryIntegrals,EI);
	end for;
	return Periods,ElementaryIntegrals;
end function;


