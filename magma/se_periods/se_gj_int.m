/*******************************************************************************

 Copyright (C) 2017 
 Code adapted from Pascal Molin,
 written by Christian Neurohr

*******************************************************************************/

// Import functions
import "se_anal_cont.m": AC_mthRoot, AC_mthRoot2, AC_mthRoot3;
import "se_help_funcs.m": MakeCCVector, PolynomialShiftVector;
import "se_spanning_tree.m": GC_AJM_Weight;
import "se_de_int.m": DE_Integrals_Factor, DE_Integrals_Factor_AJM;

// Low precision complex field
C<I> := ComplexField(20);
RPI := Real(Pi(C));

////////////////////////////////////////////////////////////
// ***** Parameters for Gauss-Jacobi integration ***** //
////////////////////////////////////////////////////////////


function GJ_Parameters(P,AB,Prec)
	D := Log(10) * Prec;
	M_r := P[1]; r := P[2];
	achr := Argcosh(r);
	if AB eq [-1/2,-1/2] then
		// Gauss-Chebyshev
		N := Ceiling((Log(10)*Prec+Log(2*RPI*M_r)+1)/(2*achr));
	else
		// Gauss-Jacobi
		I1 := 2^(&+AB+1)*Gamma(AB[1]+1)*Gamma(AB[2]+1)/Gamma(&+AB+2);
		c1 := 4*I1*M_r;
		if AB[1] eq AB[2] then
			epow := -2;
		else
			epow := -1;
		end if;
		N := Ceiling((Log(c1)+D-Log(1-Exp(achr)^epow))/(2*achr));
	end if;
	return N;
end function;

function DistanceToEllipse(x,r)
	C<I> := Parent(x);
	assert r gt 1;
	b := Sqrt(r^2-1);
	function s(t)
		return Cos(t)*Sin(t)-r*Re(x)*Sin(t)+b*Im(x)*Cos(t);
	end function;
	function sp(t)
		return (Cos(t)^2 - Sin(t)^2) - r*Re(x)*Cos(t) - b*Im(x)*Sin(t);
	end function;
	t := Re(Arccos(x));
	repeat
		nt := t;
		t -:= s(t)/sp(t);
	until Abs(t-nt) lt 10^-4;
	xr := r*Cos(t)+I*b*Sin(t);
	return xr,Abs(x-xr);
end function;
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
// ***** Gauss-Chebyshev integration ***** //
/////////////////////////////////////////////


function GC_Intergrals_Factor( VectorIntegral, Edge, SEC )
// Multiply vector integral by constants
	C<I> := SEC`ComplexField;
	oh := C!1/2;
	n := SEC`Degree[2];

	// Polynomial shift by (b+a)/(b-a)
	PolynomialShiftVector(~VectorIntegral,Edge`Data[n],SEC`Genus,0);
	
	// c1 = 2 / ((b+a)/2)^(n/2)
	n2 := Floor(n/2);
	c1 := Edge`Data[n-1]^n2;
	if n mod 2 eq 1 then
		c1 *:= Sqrt(Edge`Data[n-1]);
	end if;
	c1 := 2/c1;	

	// * -I
	if Edge`up mod 2 eq 0 then
		c1 *:= -I;
	end if;

	VectorIntegral := [ c1 * v : v in VectorIntegral ];
	// c2 = ((b+a)/2)^(i+1)
	c2 := 1;
	for k in [1..SEC`Genus] do
		c2 *:= Edge`Data[n-1];
		VectorIntegral[k] *:= c2;
	end for;
	ElemInts := [ oh * v : v in VectorIntegral ];
	return VectorIntegral, ElemInts;
end function;


function GC_Integrals(Edge,SEC)
// Gauss-Chebychev integration
	C_0 := Zero(SEC`ComplexField); 
	GJInt := SEC`IntegrationSchemes[Edge`IntSch];
	N := GJInt[1]`NPoints;
	VectorIntegral := [ C_0 : j in [1..SEC`Genus] ];
	
	if N mod 2 eq 1 then
		// Treat zero first
		N2 := Floor((N-1)/2);
		// Compute 1/y(0)
		y := 1/AC_mthRoot(0,Edge,SEC`Zetas,2,SEC`Degree[2]-2);
		VectorIntegral[1] +:= y;
	else
		N2 := Floor(N/2);
	end if;
	for k in [1..N2] do
		// Compute 1/y(x),1/y(-x)
		x := GJInt[1]`Abscissas[k];
		mx := -x;
		y1 := 1/AC_mthRoot(x,Edge,SEC`Zetas,2,SEC`Degree[2]-2);
		y2 := 1/AC_mthRoot(mx,Edge,SEC`Zetas,2,SEC`Degree[2]-2);
		// Compute all differentials
		y1xi := y1;
		y2xi := y2;
		VectorIntegral[1] +:= y1xi + y2xi;
		for j in [2..SEC`Genus] do
			y1xi *:= x;
			y2xi *:= mx;
			VectorIntegral[j] +:= y1xi + y2xi;
		end for;
	end for;

	// Multiply by weight Pi/n
	for j in [1..SEC`Genus] do
		VectorIntegral[j] *:= GJInt[1]`Weights[1];
	end for;
	return VectorIntegral;
end function;


function GC_Integrals_Edge(Edge,SEC)
// Integrate an edge of the spanning tree with Gauss-Chebychev integration

	// Compute integral on [-1,1]
	VectorIntegral := GC_Integrals(Edge,SEC);

	// Multiplication with constants
	PeriodsEdge, ElemIntegralEdge := GC_Intergrals_Factor(VectorIntegral,Edge,SEC);

	return PeriodsEdge, ElemIntegralEdge; 
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

	assert &and[ AB[1] gt -1, AB[2] gt -1];

	// Careful: AB[1] corresponds to (1-u) and AB[2] to (1+u)
	Reverse(~AB);
	Abscissas := []; 
	Weights := [];
	R := RealField(MaxPrec);

	if AB[1] eq AB[2] then
		n := Floor((N+1)/2);
	else
		n := N;
	end if;

	// Gauss-Chebyshev
	if AB eq [-1/2,-1/2] then
		const := 1/(2*N)*Pi(R);
		for k in [1..n] do
			Append(~Abscissas,Cos((2*k-1)*const));
		end for;
		Append(~Weights,2*const);
		return Abscissas,Weights;
	end if;

	// From here on: Gauss-Jacobi

	// Error
	Eps := Real(10^-(MaxPrec-4));
	
	// Precomputations
	RN := R!N; N2 := N^2;	
	B := []; C := [];
	APB := R!(AB[1]+AB[2]);
	AMB := R!(AB[1]-AB[2]);
	NAB := R!(2*(N+AB[1])*(N+AB[2]));
	ASQMBSQ := AB[1]^2 - AB[2]^2;
	for j in [1..N] do
		tmp := 2*j + APB;
		tmp2 := 1/(2*j*(j+APB)*(tmp-2));
		Append(~B,[tmp*(tmp-1)*(tmp-2)*tmp2,(tmp-1)*ASQMBSQ*tmp2]);
		Append(~C,2*(j-1+AB[1])*(j-1+AB[2])*tmp*tmp2);
	end for;
	//ChangeUniverse(~A,R);
	OH := R!(1/2);
	AMBOT := AMB * OH;
	APBP2OT := (APB+2) * OH;
	NAMB := N*AMB;
	Ntmp := N*tmp;
	Const := tmp*2^APB*Gamma(AB[1]+RN)*Gamma(AB[2]+RN)/(Gamma(RN+1)*Gamma(RN+APB+1));

	// Start with low precision
	sp := 16;
	RL := RealField(sp);
	
	for k in [1..n] do
		if k eq 1 then
			an := AB[1]/N;
			bn := AB[2]/N;
			r1 := (1+AB[1])*((139/50)/(4+N2) + (98/125)*an/N);
			r2 := 1+(37/25)*an+(24/25)*bn+(113/250)*an^2+(83/100)*an*bn;
			z := 1-r1/r2;
		elif k eq 2 then
			r1 := (41/10+AB[1])/((1+AB[1])*(1+(39/250)*AB[1]));
			r2 := 1+3/50*(N-8)*(1+3/25*AB[1])/N;
			r3 := 1+3/250*AB[2]*(1+(1/4)*Abs(AB[1]))/N;
			z -:= (1-z)*r1*r2*r3;
		elif k eq 3 then
			r1 := (1.67+7/25*AB[1])/(1+37/100*AB[1]);
			r2 := 1+11/50*(N-8)/N;
			r3 := 1+8*AB[2]/((157/25+AB[2])*N2);
			z -:= (Abscissas[1]-z)*r1*r2*r3;
		elif k eq N-1 then
			r1 := (1+47/200*AB[2])/(383/500+119/1000*AB[2]);
			r2 := 1/(1+639/1000*(N-4)/(1+71/100*(N-4)));
			r3 := 1/(1+20*AB[1]/((15/2+AB[1])*N2));
			z +:= (z-Abscissas[N-3])*r1*r2*r3;
		elif k eq N then
			r1 := (1+37/100*AB[2])/(167/100+7/25*AB[2]);
			r2 := 1/(1+11/50*(N-8)/N);
			r3 := 1/(1+8*AB[1]/((157/25+AB[1])*N2));
			z +:= (z-Abscissas[N-2])*r1*r2*r3;
		else
			z := 3*Abscissas[k-1]-3*Abscissas[k-2]+Abscissas[k-3];
		end if;
		z := RL!z;
		p := sp;
    		repeat
			p := Min(2*p,MaxPrec);
			z := ChangePrecision(z,p);
	     		//p1,ppinv := RS_JacobiPolyEval(z,AB,N);
			//tmp := 2+APB;
			p1 := AMBOT+APBP2OT*z;
			p2 := 1;
			for j in [2..N] do
				p3 := p2;
				p2 := p1;
				//tmp := 2*j + APB;
				//a := 2*j*(j+APB)*(tmp-2);
				//b := (tmp-1)*(AB[1]^2-AB[2]^2+tmp*(tmp-2)*z);
				//b := B[j][1]*z + B[j][2];
				//c := 2*(j-1+AB[1])*(j-1+AB[2])*tmp;
				//p1 := (b*p2-c*p3)/a;
				//p1 := ((B[j][1]*z + B[j][2])*p2-C[j]*p3)*A[j];
				p1 := (B[j][1]*z + B[j][2])*p2-C[j]*p3;
			end for;
			ppinv := tmp*(1-z*z)/((NAMB-Ntmp*z)*p1+NAB*p2);
      			fxk := p1*ppinv; 
			z -:= fxk;
		until Abs(fxk) lt Eps;
		Abscissas[k] := z;
		//Weights[k] := Const*tmp*ppinv/p2;
		Weights[k] := Const*ppinv/p2;
    	end for;
	return Abscissas,Weights;
end function;


// Define type GJ_Int
declare type GJ_Int;
declare attributes GJ_Int: Abscissas, Weights, AlphaBeta, NPoints, Degree, Params, Prec;


// Constructor
intrinsic GJ_Integration( Params, SEC : AJM := false )
{ Construct the double exponential integration scheme }

	// Parameters = < M_i, r_i >
	// N_i < N_{i+1}
	m := SEC`Degree[1];
	DFF := SEC`HolomorphicDifferentials;
	ALLGJ_Integrations := [];
	// For each N
	if AJM then
		NSchemes := #SEC`IntegrationSchemes; 
		if NSchemes eq 0 or SEC`IntegrationSchemes[1][1]`Prec lt Precision(SEC`ComplexField) then 
			for P in Params do
				// Compute m-1 Schemes
				GJ_Integrations := [];
				N := GJ_Parameters(P,[-(m-1)/m,0/1],SEC`Prec);
				for j in [DFF[1]..m-1] do
					GJInt := New(GJ_Int);
					GJInt`Degree := m;
					GJInt`AlphaBeta := [-j/m,0/1];
					GJInt`NPoints := N;
					GJInt`Params := P;
					GJInt`Abscissas, GJInt`Weights := GJ_IntegrationPoints(GJInt`AlphaBeta,GJInt`NPoints,Precision(SEC`ComplexField));
					GJInt`Prec := Precision(SEC`ComplexField);
					Append(~GJ_Integrations,GJInt);
				end for;
				Append(~ALLGJ_Integrations,GJ_Integrations);
			end for;
			SEC`IntegrationSchemes := ALLGJ_Integrations;
		else
			Ns := [ IntSch[1]`NPoints : IntSch in SEC`IntegrationSchemes ];
			for k in [1..#Params] do
				// Compute m-1 Schemes
				GJ_Integrations := [];
				N := GJ_Parameters(Params[k],[-(m-1)/m,0/1],SEC`Prec);
				if N gt Ns[1] then
					for j in [DFF[1]..m-1] do
						GJInt := New(GJ_Int);
						GJInt`Degree := m;
						GJInt`AlphaBeta := [-j/m,0/1];
						GJInt`NPoints := N;
						GJInt`Params := Params[k];
						GJInt`Abscissas, GJInt`Weights := GJ_IntegrationPoints(GJInt`AlphaBeta,GJInt`NPoints,Precision(SEC`ComplexField));
						GJInt`Prec := Precision(SEC`ComplexField);
						Append(~GJ_Integrations,GJInt);
					end for;
					Append(~ALLGJ_Integrations,GJ_Integrations);
					Ns[1] := N;
				else
					l := Max( [ j : j in [1..NSchemes] | Ns[j] ge N ] );
					Append(~ALLGJ_Integrations,SEC`IntegrationSchemes[l]);
				end if;
			end for;
		end if;
	else
		for P in Params do
			// m-1 Schemes
			GJ_Integrations := [];
			N := GJ_Parameters(P,[-(m-1)/m,-(m-1)/m],SEC`Prec);
			for j in [DFF[1]..m-1] do
				GJInt := New(GJ_Int);
				GJInt`Degree := m;
				GJInt`AlphaBeta := [-j/m,-j/m];
				GJInt`NPoints := N;
				GJInt`Params := P;
				GJInt`Abscissas, GJInt`Weights := GJ_IntegrationPoints(GJInt`AlphaBeta,GJInt`NPoints,Precision(SEC`ComplexField));
				if GJInt`AlphaBeta eq [-1/2,-1/2] and m gt 2 then
					GJInt`Weights := [ GJInt`Weights[1] : k in [1..Floor(N/2)+1] ];
				end if; 
				GJInt`Prec := Precision(SEC`ComplexField);
				Append(~GJ_Integrations,GJInt);
			end for;
			Append(~ALLGJ_Integrations,GJ_Integrations);
		end for;
		SEC`IntegrationSchemes := ALLGJ_Integrations;
	end if;
	
end intrinsic;


// Printing
intrinsic Print(GJInt::GJ_Int)
{ Printing }
	print "Params:",GJInt`Params;
	print "NPoints:",GJInt`NPoints;
	print "[Alpha,Beta]:",GJInt`AlphaBeta;
	print "Degree:",GJInt`Degree;
	print "Precision:",GJInt`Prec;
end intrinsic;


// Gauss-Jacobi integration

function GJ_Integrals(Edge,SEC)
	C_0 := Zero(SEC`ComplexField); 
	GJInts := SEC`IntegrationSchemes[Edge`IntSch];
	N := GJInts[1]`NPoints;
	DFF := SEC`HolomorphicDifferentials;
	VectorIntegral := [ C_0 : j in [1..SEC`Genus] ];
	N2 := Floor(N/2);
	
	if N mod 2 eq 1 then
		ct := 1;
		y := 1/AC_mthRoot(0,Edge,SEC`Zetas,SEC`Degree[1],SEC`Degree[2]-2); // 1/y(0)
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
			mx := -x;
			y1 := 1/AC_mthRoot(x,Edge,SEC`Zetas,SEC`Degree[1],SEC`Degree[2]-2); // 1/y(x)
			y2 := 1/AC_mthRoot(mx,Edge,SEC`Zetas,SEC`Degree[1],SEC`Degree[2]-2); // 1/y(-x)
			wy1 := GJInts[j]`Weights[t] * y1^(DFF[4][ct]);
			wy2 := GJInts[j]`Weights[t] * y2^(DFF[4][ct]);
			VectorIntegral[ct] +:= wy1 + wy2;
			ct +:= 1;
			for k in [1..DFF[3][j]-1] do
				wy1 *:= x;
				wy2 *:= mx;
				VectorIntegral[ct] +:= wy1 + wy2;
				ct +:= 1;
			end for;
		end for;
	end for;
	return VectorIntegral;
end function;


function GJ_Integrals_Edge(EdgeData,SEC)
// Integrate an edge of the spanning tree with Gauss-Chebychev integration

	// Compute integral on [-1,1]
	VectorIntegral := GJ_Integrals(EdgeData,SEC);

	// Multiplication with constants
	PeriodsEdge, ElemIntegralEdge := DE_Integrals_Factor(VectorIntegral,EdgeData,SEC);

	return PeriodsEdge, ElemIntegralEdge; 
end function;


function GJ_Integrals_Tree(SEC)
// Compute integrals for spanning tree
	Periods := []; 
	ElementaryIntegrals := [];
	for k in [1..SEC`SpanningTree`Length] do
		vprint SE,2 : "Integrating edge #:",k;
		if SEC`Degree[1] eq 2 then
			Period, EI := GC_Integrals_Edge(SEC`SpanningTree`Edges[k],SEC);
			Append(~Periods,[ [ P ] : P in Period ]);
		else
			Period, EI := GJ_Integrals_Edge(SEC`SpanningTree`Edges[k],SEC);
			Append(~Periods,[ P : P in Period ]);
		end if;
		Append(~ElementaryIntegrals,EI);
	end for;
	return Periods,ElementaryIntegrals;
end function;


//////////////////////////////////////////////////////////
// ***** Numerical Integration of Abel-Jacobi map ***** //
//////////////////////////////////////////////////////////

function GJ_Integrals_AJM(Edge,SEC)
	C_0 := Zero(SEC`ComplexField); 
	GJInts := SEC`IntegrationSchemes[Edge`IntSch];
	DFF := SEC`HolomorphicDifferentials;
	VectorIntegral := [ C_0 : j in [1..SEC`Genus] ];

	// Evaluate differentials at abisccsas
	for t in [1..GJInts[1]`NPoints] do
		ct := 1;
		for j in [1..DFF[2]] do
			x := GJInts[j]`Abscissas[t];
			y := 1/AC_mthRoot(x,Edge,SEC`Zetas,SEC`Degree[1],SEC`Degree[2]-1); // 1/y(x)
			wy := GJInts[j]`Weights[t] * y^(DFF[4][ct]);
			VectorIntegral[ct] +:= wy;
			ct +:= 1;
			for k in [1..DFF[3][j]-1] do
				wy *:= x;
				VectorIntegral[ct] +:= wy;
				ct +:= 1;
			end for;
		end for;
	end for;
	return VectorIntegral;
end function;

function GJ_Integrals_Edge_AJM(Edge,SEC)
// Integrate an "complex edge" for the Abel-Jacobi map

	// Numerical integration on [-1,1]
	VectorIntegral := GJ_Integrals_AJM(Edge,SEC);

	// Multiplication with constants (same as for DE)
	Integral_Edge_AJM := DE_Integrals_Factor_AJM(VectorIntegral,Edge,SEC);

	return Integral_Edge_AJM;
end function;
