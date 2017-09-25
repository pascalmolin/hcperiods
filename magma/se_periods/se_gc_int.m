/*******************************************************************************

 Copyright (C) 2017 
 Code adapted from Pascal Molin,
 written by Christian Neurohr

*******************************************************************************/

// Import functions
import "se_anal_cont.m": AC_mthRoot;
import "se_help_funcs.m": MakeCCVector, PolynomialShiftVector;
import "se_spanning_tree.m": GC_AJM_Weight;
import "se_de_int.m": DE_Integrals_Factor;

// Low precision complex field
C<I> := ComplexField(40);
RPI := Real(Pi(C));

////////////////////////////////////////////////////////////
// ***** Parameters for Gauss-Chebychev integration ***** //
////////////////////////////////////////////////////////////

/*function UpperBound_M(M_0,r,V_r,Len)
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
end function;*/



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
	if N mod 2 eq 1 then
		// Treat zero first
		N2 := Floor((N-1)/2);
		// Compute 1/y(0)
		y := 1/AC_mthRoot(0,EdgeData,SEC`Zetas,up,2,n-2);
		VectorIntegral[1] +:= y;
	else
		N2 := Floor(N/2);
	end if;
	for k in [1..N2] do
		// Compute x
		x := Cos((2*k-1)*const);

		// Compute 1/y(x)
		y1 := 1/AC_mthRoot(x,EdgeData,SEC`Zetas,up,2,n-2);
		y2 := 1/AC_mthRoot(-x,EdgeData,SEC`Zetas,up,2,n-2);

		// Compute all differentials
		y1xi := y1;
		y2xi := y2;
		VectorIntegral[1] +:= y1xi + y2xi;
		for j in [2..g] do
			y1xi *:= x;
			y2xi *:= -x;
			VectorIntegral[j] +:= y1xi + y2xi;
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
	Params := SEC`SpanningTree`Params;
	Periods := []; ElementaryIntegrals := [];
	for k in [1..SEC`SpanningTree`Length] do
		Period, EI := GC_Integrals_Edge(SEC`SpanningTree`Data[k],SEC,Params[Round(SEC`SpanningTree`Edges[k][3])][3]);
		Append(~Periods,[ [ P ] : P in Period ]);
		Append(~ElementaryIntegrals,EI);
	end for;
	return Periods,ElementaryIntegrals;
end function;




/////////////////////////////////////////////
// ***** Gauss-Jacobi integration ***** //
/////////////////////////////////////////////


function JacobiPolyEval(z,AB,N)
// Evaluate the N-th Legendre polynomial and its derivative at x using the three-term recurrence relation
	APB := AB[1]+AB[2];
	tmp := 2+APB;
	p1 := (AB[1]-AB[2]+tmp*z)/2;
	p2 := 1;
	for j in [2..N] do
		p3 := p2;
		p2 := p1;
		tmp := 2*j + APB;
		a := 2*j*(j+APB)*(tmp-2);
		b := (tmp-1)*(AB[1]^2-AB[2]^2+tmp*(tmp-2)*z);
		c := 2*(j-1+AB[1])*(j-1+AB[2])*tmp;
		p1 := (b*p2-c*p3)/a;
	end for;
	//pp := (N*(AB[1]-AB[2]-tmp*z)*p1+2*(N+AB[1])*(N+AB[2])*p2)/(tmp*(1-z^2));
	ppinv := (tmp*(1-z^2))/(N*(AB[1]-AB[2]-tmp*z)*p1+2*(N+AB[1])*(N+AB[2])*p2);
	return p1,ppinv;
end function;

function GJ_IntegrationPoints( AB, N, MaxPrec )
// Computes Gauss-Jacobi quadrature on N points of weight function AB[1],AB[2] to precision MaxPrec-5
	R := RealField(MaxPrec);
	RN := R!N; N2 := N^2;
	
	Eps := Real(10^-(MaxPrec-4));
	alpha := AB[1];
	beta := AB[2];
	APB := AB[1]+AB[2];
	Const := R!2^APB*Gamma(alpha+RN)*Gamma(beta+RN)/(Gamma(RN+1)*Gamma(RN+alpha+beta+1));
	assert &and[ alpha gt -1, beta gt -1];

	Abscissas := []; Weights := [];
	sp := 16;
	RL := RealField(sp);
	
	if alpha eq beta then
		n := Floor((N+1)/2);
	else
		n := N;
	end if;

	for k in [1..n] do
		if k eq 1 then
			an := alpha/N;
			bn := beta/N;
			r1 := (1+alpha)*((139/50)/(4+N2) + (98/125)*an/N);
			r2 := 1+(37/25)*an+(24/25)*bn+(113/250)*an^2+(83/100)*an*bn;
			z := 1-r1/r2;
		elif k eq 2 then
			r1 := (41/10+alpha)/((1+alpha)*(1+(39/250)*alpha));
			r2 := 1+3/50*(N-8)*(1+3/25*alpha)/N;
			r3 := 1+3/250*beta*(1+(1/4)*Abs(alpha))/N;
			z -:= (1-z)*r1*r2*r3;
		elif k eq 3 then
			r1 := (1.67+7/25*alpha)/(1+37/100*alpha);
			r2 := 1+11/50*(N-8)/N;
			r3 := 1+8*beta/((157/25+beta)*N2);
			z -:= (Abscissas[1]-z)*r1*r2*r3;
		elif k eq N-1 then
			r1 := (1+47/200*beta)/(383/500+119/1000*beta);
			r2 := 1/(1+639/1000*(N-4)/(1+71/100*(N-4)));
			r3 := 1/(1+20*alpha/((15/2+alpha)*N2));
			z +:= (z-Abscissas[N-3])*r1*r2*r3;
		elif k eq N then
			r1 := (1+37/100*beta)/(167/100+7/25*beta);
			r2 := 1/(1+11/50*(N-8)/N);
			r3 := 1/(1+8*alpha/((157/25+alpha)*N2));
			z +:= (z-Abscissas[N-2])*r1*r2*r3;
		else
			z := 3*Abscissas[k-1]-3*Abscissas[k-2]+Abscissas[k-3];
		end if;
		z := RL!z;
		alphabeta := alpha+beta;
		p := sp;
    		repeat
			p := Min(2*p,MaxPrec);
			z := ChangePrecision(z,p);
	     		//p1,ppinv := RS_JacobiPolyEval(z,AB,N);
			tmp := 2+APB;
			p1 := (AB[1]-AB[2]+tmp*z)/2;
			p2 := 1;
			for j in [2..N] do
				p3 := p2;
				p2 := p1;
				tmp := 2*j + APB;
				a := 2*j*(j+APB)*(tmp-2);
				b := (tmp-1)*(AB[1]^2-AB[2]^2+tmp*(tmp-2)*z);
				c := 2*(j-1+AB[1])*(j-1+AB[2])*tmp;
				p1 := (b*p2-c*p3)/a;
			end for;
			ppinv := (tmp*(1-z^2))/(N*(AB[1]-AB[2]-tmp*z)*p1+2*(N+AB[1])*(N+AB[2])*p2);
      			nz := z;
			z -:= p1*ppinv;
		until Abs(z-nz) lt Eps;
		Abscissas[k] := z;
		Weights[k] := Const*tmp*ppinv/p2;
    	end for;

	if alpha eq beta then
		return Abscissas,Weights;
	else
		return Abscissas,Weights;
	end if;
end function;

// Define type GJ_Int
declare type GJ_Int;
declare attributes GJ_Int: Abscissas, Weights, AlphaBeta, NPoints, Degree, Params;


// Constructor
intrinsic GJ_Integration( Params, SEC ) -> SeqEnum[SeqEnum[GJ_Int]]
{ Construct the double exponential integration scheme }

	// Parameters = < M_i, r_i, N_i >
	// N_i < N_{i+1}
	m := SEC`Degree[1];
	DFF := SEC`HolomorphicDifferentials;
	ALLGJ_Integrations := [];
	// For each N
	for P in Params do
		// m-1 Schemes
		GJ_Integrations := [];
		for j in [DFF[1]..m-1] do
			GJInt := New(GJ_Int);
			GJInt`Degree := m;
			GJInt`Params := P;
			GJInt`AlphaBeta := [-j/m,-j/m];
			GJInt`Abscissas, GJInt`Weights := GJ_IntegrationPoints(GJInt`AlphaBeta,GJInt`Params[3],Precision(SEC`ComplexField));
			Append(~GJ_Integrations,GJInt);
		end for;
		Append(~ALLGJ_Integrations,GJ_Integrations);
	end for;
	return ALLGJ_Integrations;
end intrinsic;


// Printing
intrinsic Print(GJInt::GJ_Int)
{ Printing }
	print "Params:",GJInt`Params;
	print "[Alpha,Beta]:",GJInt`AlphaBeta;
	print "Degree:",GJInt`Degree;
end intrinsic;


// Gauss-Jacobi integration

function GJ_Integrals(EdgeData,SEC,GJInts)
	C<I> := SEC`ComplexField;
	C_0 := Zero(C); g := SEC`Genus; m := SEC`Degree[1]; n := SEC`Degree[2]; N := GJInts[1]`Params[3];
	DFF := SEC`HolomorphicDifferentials;
	VectorIntegral := [ C_0 : j in [1..g] ];
	up := Floor(EdgeData[n+1]);
	N2 := Floor(N/2);

	if N mod 2 eq 1 then
		ct := 1;
		y := 1/AC_mthRoot(0,EdgeData,SEC`Zetas,up,m,n-2); // 1/y(0)
		for j in [1..DFF[2]] do
			wy := GJInts[j]`Weights[N2+1] * y^(DFF[4][ct]);
			VectorIntegral[ct] +:= wy;
			ct +:= DFF[3][j];
		end for;
	end if;

	// Evaluate differentials at abisccsas

	for t in [1..N2] do
		ct := 1;
		for j in [1..DFF[2]] do
			x := GJInts[j]`Abscissas[t];
			y1 := 1/AC_mthRoot(x,EdgeData,SEC`Zetas,up,m,n-2); // 1/y(x)
			y2 := 1/AC_mthRoot(-x,EdgeData,SEC`Zetas,up,m,n-2); // 1/y(-x)
			wy1 := GJInts[j]`Weights[t] * y1^(DFF[4][ct]);
			wy2 := GJInts[j]`Weights[t] * y2^(DFF[4][ct]);
			VectorIntegral[ct] +:= wy1 + wy2;
			ct +:= 1;
			for k in [1..DFF[3][j]-1] do
				wy1 *:= x;
				wy2 *:= -x;
				VectorIntegral[ct] +:= wy1 + wy2;
				ct +:= 1;
			end for;
		end for;
	end for;
	return VectorIntegral;
end function;



function GJ_Integrals_Edge( EdgeData, SEC, GJInts)
// Integrate an edge of the spanning tree with Gauss-Chebychev integration

	// Compute integral on [-1,1]
	VectorIntegral := GJ_Integrals(EdgeData,SEC,GJInts);

	// Multiplication with constants
	PeriodsEdge, ElemIntegralEdge := DE_Integrals_Factor(VectorIntegral,EdgeData,SEC);

	return PeriodsEdge, ElemIntegralEdge; 
end function;


function GJ_Integrals_Tree(SEC, ALLGJInts)
// Compute integrals for tree
	// Compute integration parameters for spanning tree
	Params := SEC`SpanningTree`Params;
	Periods := []; ElementaryIntegrals := [];
	for k in [1..SEC`SpanningTree`Length] do
		P, EI := GJ_Integrals_Edge(SEC`SpanningTree`Data[k],SEC,ALLGJInts[Round(SEC`SpanningTree`Edges[k][3])]);
		Append(~Periods,P);
		Append(~ElementaryIntegrals,EI);
	end for;
	return Periods,ElementaryIntegrals;
end function;
