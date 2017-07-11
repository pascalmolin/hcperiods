/*******************************************************************************

 Copyright (C) 2017 
 Code adapted from Pascal Molin,
 written by Christian Neurohr

*******************************************************************************/

// Import functions
import "se_anal_cont.m": AC_mthRoot;
import "se_help_funcs.m": MakeCCVector, PolynomialShiftVector;
import "se_spanning_tree.m": GC_AJM_Weight;

// Low precision complex field
C30<I> := ComplexField(30);
Pi30 := Real(Pi(C30));

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
	ChangeUniverse(~CCV,C30);
	for k in [1..Len] do
		r_k := (1/2) * ( Abs(CCV[k]+1) + Abs(CCV[k]-1) );
		assert r_k gt 1;
		r_0 := Min(r_0,r_k);
		Append(~V_r,r_k);
	end for;
	// M_0 = D + Log(2Pi) + n*Log(r_0)
	M_0 := Log(2*Pi30) + n*Log(r_0) + Log(10)*Prec;
	vprint SE,2 : "M_0:",M_0;

   	// Make it precise
	r := (1/2) * (r_0 + 1);

	// r small enough?
	K,M_r := UpperBound_M(M_0,r,V_r,Len);
	while K ne 0 do
		A := 1 + 2*M_r/K;
		r := Cosh(r_0 - r_0/(Log(A/Argcosh(r_0))+A));
		K,M_r := UpperBound_M(M_0,r,V_r,Len);
	end while;


	vprint SE,2 : "r:",r;
	vprint SE,2 : "M(r):",M_r;

	// Number of integration points
	r := Argcosh(r); // r' = Log(r+Sqrt(r^2-1))
	M_r := (M_r/(2*r));
	N := Ceiling(M_r);
	//N := Ceiling(M_r/(2*Argcosh(r)));

	vprint SE,2 : "N:",N;

	// Error
	//e :=  Exp(M_r)/(Exp(2*N*r)-1);

	return N;
end function;
function GC_Params_Tree( STree,Prec);
	return GC_Params(STree`Data[STree`WorstEdge],STree`Length-1,STree`Length+1,Prec);
end function;


/////////////////////////////////////////////
// ***** Gauss-Chebychev integration ***** //
/////////////////////////////////////////////


function GC_Intergrals_Factor( VectorIntegral, EdgeData, SEC )
// Multiply vector integral by constants
	
	C<I> := SEC`ComplexField;
	g := SEC`Genus;
	n := SEC`Degree[2];
	up := Floor(EdgeData[n+1]);

	// Polynomial shift by (b+a)/(b-a)
	PolynomialShiftVector( ~VectorIntegral, EdgeData[n], g, 0 );
	
	// c1 = 2 / ((b+a)/2)^(n/2)
	n2 := Floor(n/2);
	c1 := EdgeData[n-1]^n2;
	if n mod 2 eq 1 then
		c1 *:= Sqrt(EdgeData[n-1]);
	end if;
	c1 := 2/c1;	

	// * -I
	if up mod 2 eq 0 then
		c1 *:= -I;
	end if;

	VectorIntegral := [ c1 * v : v in VectorIntegral ];
	// c2 = ((b+a)/2)^(i+1)
	c2 := 1;
	for k in [1..g] do
		c2 *:= EdgeData[n-1];
		VectorIntegral[k] *:= c2;
	end for;
	ElemInts := [ (1/2) * v : v in VectorIntegral ];
	return VectorIntegral, ElemInts;
end function;


// EdgeData = [ u_1, ... , u_{n-2}, (b-a)/2, (b+a)/(b-a), up ]
function GC_Integrals( EdgeData, SEC, N )
// Gauss-Chebychev integration
	C<I> := SEC`ComplexField;
	g := SEC`Genus;
	n := SEC`Degree[2];
	VectorIntegral := [ Zero(C) : j in [1..g] ]; 
	PiR := Real(Pi(C));
	const := 1/(2*N)*PiR;
	up := Floor(EdgeData[n+1]);
	for k in [1..N] do
		// Compute x
		x := Cos((2*k-1)*const);

		// Compute 1/y(x)
		y := 1/AC_mthRoot(x,EdgeData,SEC`Zetas,up,2,n-2);

		// Compute all differentials
		yxi := y;
		VectorIntegral[1] +:= yxi;
		for j in [2..g] do
			yxi *:= x;
			VectorIntegral[j] +:= yxi;
		end for;
	end for;

	// Multiply by weight Pi/n
	w := PiR/N;
	for j in [1..g] do
		VectorIntegral[j] *:= w;
	end for;
	return VectorIntegral;
end function;


function GC_Integrals_Edge( EdgeData, SEC, N)
// Integrate an edge of the spanning tree with Gauss-Chebychev integration

	// Compute integral on [-1,1]
	VectorIntegral := GC_Integrals(EdgeData,SEC,N);

	// Multiplication with constants
	PeriodsEdge, ElemIntegralEdge := GC_Intergrals_Factor(VectorIntegral,EdgeData,SEC);

	return PeriodsEdge, ElemIntegralEdge; 
end function;


function GC_Integrals_Tree(SEC)
// Compute integrals for tree
	// Compute integration parameters for spanning tree
	NPoints := SEC`SpanningTree`Params[2]; NInts := #NPoints;

	Periods := []; 
	ElementaryIntegrals := [];
	for k in [1..SEC`SpanningTree`Length] do
		l := NInts;
		while SEC`SpanningTree`Edges[k][3] lt SEC`SpanningTree`Params[3][l] do
			l -:= 1;
		end while;
		Period, EI := GC_Integrals_Edge(SEC`SpanningTree`Data[k],SEC,NPoints[l]);
		Append(~Periods,[ [ P ] : P in Period ]);
		Append(~ElementaryIntegrals,EI);
	end for;
	return Periods,ElementaryIntegrals;
end function;
