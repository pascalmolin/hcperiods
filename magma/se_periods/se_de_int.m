/******************************************************************************

 Copyright (C) 2017 
 Code adapted from Pascal Molin,
 written by Christian Neurohr

 ******************************************************************************/

import "se_anal_cont.m": AC_mthRoot, AC_mthRoot2, AC_mthRoot3;
import "se_help_funcs.m": PolynomialShiftVector;
import "se_spanning_tree.m": DE_AJM_Weight;


///////////////////////////////////////////////////////////////
// ***** Parameters for double-exponential integration ***** //
///////////////////////////////////////////////////////////////


C<I> := ComplexField(20);
RPI := Real(Pi(C));

function Distance_1(P)
	xP := Abs(Re(P));
	if xP gt 1 then
		return Abs(xP-1+Parent(P).1*Im(P));
	else
		return Abs(Im(P));
	end if;
end function;

function Bound_M1(CCV,len,m : AJM := false)
// Compute M1
	M1 := (1/&*[ Distance_1(CCV[k]) : k in [1..len] ])^(1/m);
	if AJM then
		M1 *:= 2^(1/m);
	end if;
	if M1 gt 1 then
		return M1^(m-1);
	else
		return M1;
	end if;
end function;


function Bound_M2(CCV,len,m,n,r : Lambda := RPI/2, AJM := false )
// Compute the bound M2 (approximation)

	// Search on interval [0,t_max]
	tmax := Argcosh(RPI/(2*Lambda*Sin(r)));
	//print "tmax:",tmax;

	Phi := function(t)
		return Tanh(Lambda*Sinh(t+I*r));
	end function;
	Dist := function(z,P)
		P := Abs(Re(P)) + I*Abs(Im(P));
		return Abs(z-P);
	end function;

	// Stepsize
	ss := 20;

	// Abs(1-tanh(sinh(I*r)))
	MDTO := Abs(1-Phi(0));

	M2 := 1;
	for k in [0..ss] do
		t := k/ss*tmax;
		z := Phi(t);
		za := Abs(z);
		if za gt 1 then
			km := za^(n-2);
		else
			km := za;
		end if;
		if AJM then
			u := (MDTO/&*[ Dist(z,CCV[l]): l in [1..len] ])^(1/m);
		else
			u := (1/&*[ Dist(z,CCV[l]): l in [1..len] ])^(1/m);
		end if;
		if u gt 1 then
			u := u^(m-1);
		end if;
		km *:= u;
		if km gt M2 then
			xm := t;
		end if;
		M2 := Max(km,M2);
	end for;
	return M2;
end function;



/////////////////////////////////////////////////////
// ***** Double-exponential integration type ***** //
/////////////////////////////////////////////////////


// Define  type DE_Int
declare type DE_Int;
declare attributes DE_Int: r, Factor, Steplength, Lambda, Abscissas, Weights, NPoints, Degree, ExtraFactors, Prec;


// Integration parameters
procedure DE_IntegrationPoints( DEInt : AJM := false )
// Computes integration points and weights for double-exponential integration on the interval [-1,1]

	h := DEInt`Steplength;
	R := Parent(h);
	mi := 1/DEInt`Degree;
	oh := R!(1/2);
	DEInt`Abscissas := [Zero(R)];
	DEInt`Weights := [<One(R),One(R)>];
	
	eh := Exp(h); // e^h
	eh_inv := 1/eh; // e^(-h)
	ekh := 1; // e^(0h)
	ekh_inv := 1; // e^(0-h)

	for k in [1..DEInt`NPoints] do
		ekh *:= eh; // e^(kh)
		ekh_inv *:= eh_inv; // e^(-kh)
     		sh := oh*(ekh-ekh_inv); // sinh(kh) = (1/2) * (e^(kh) - e^(-kh))
		ch := sh+ekh_inv; ; // cosh(kh) = (e^(kh) + e^(-kh))
      		esh := Exp(DEInt`Lambda*sh); // e^(Lambda*sinh(kh))
      		esh_inv := 1/esh; // e^-(Lambda*sinh(kh))
		chsh := oh*(esh+esh_inv); // cosh(Lambda*sinh(kh)))
		shsh := chsh-esh_inv; // sinh(Lambda*sinh(kh))
		Append(~DEInt`Abscissas,shsh/chsh); // tanh(Lambda*sinh(kh)) =  sinh(Lambda*sinh(kh)) / cosh(Lambda*sinh(kh))
		chsh *:= chsh; // cosh(Lambda*sinh(kh))^2
		// <cosh(kh)/cosh(Lambda*sinh(kh))^2,cosh(Lambda*sinh(kh))^(2/m)>
		Append(~DEInt`Weights,<ch/chsh,chsh^mi>);
     	end for;
	if AJM then
		DEInt`ExtraFactors := [<One(R),One(R)>] cat [ <(1-thsh)^mi,(1+thsh)^mi> : thsh in Remove(DEInt`Abscissas,1) ];
	end if;
end procedure;

// Constructor
intrinsic DE_Integration( Params, SEC : AJM := false )
{ Construct the double exponential integration scheme }

	// Parameters = < M1,M2,[ r_min,..,r_max ]>
	// r_i < r_{i+1}

	R := RealField(Precision(SEC`ComplexField));
	RPi2 := Pi(R)/2;
	Lambda := RPi2;
	m := SEC`Degree[1];

	// 'Worst' Alpha
	Alpha := 1/m;

	// Compute D
	D := SEC`Prec * Log(R!10);

	SEC`IntegrationSchemes := [];
	for P in Params do
		// Get parameters;
		M1 := P[1];M2 := P[2];r := R!P[3];
		// Compute N,h
		Alpha := 1/m;
		X_r := Cos(r) * Sqrt( ( Pi(R) / (2*Lambda*Sin(r))) - 1 );
		B_ra := (2/Cos(r)) * ( (X_r/2) * (1/(Cos(Lambda*Sin(r))^(2*Alpha)) + (1/X_r^(2*Alpha))) ) + (1/(2*Alpha*Sinh(X_r)^(2*Alpha)));
		h := Real( 4 * RPi2 * r /  ( D+Log(2*M2 * B_ra + 1)));
		N := Ceiling(Argsinh((D+ Log((2^(2*Alpha+1)*M1)/Alpha )) / ( 2*Alpha*Lambda ))/h);
		// New scheme
		DEInt := New(DE_Int);
		DEInt`Prec := Precision(R); 
		DEInt`r := ChangePrecision(r,10);
		DEInt`Steplength := h;
		DEInt`NPoints := N;
		DEInt`Lambda := Lambda;
		DEInt`Degree := m;
		DEInt`Factor := DEInt`Lambda * DEInt`Steplength;	
		// Compute weights and abscissas
		DE_IntegrationPoints(DEInt:AJM:=AJM);
		Append(~SEC`IntegrationSchemes,DEInt);
	end for;
end intrinsic;


// Printing
intrinsic Print(DEInt::DE_Int)
{ print }
	print "r:",ChangePrecision(DEInt`r,10);
	print "Steplength:",ChangePrecision(DEInt`Steplength,10);
	print "Number of abscissas:",DEInt`NPoints;
	print "Lambda:",ChangePrecision(DEInt`Lambda,10);
	print "Degree:",DEInt`Degree;
	print "Factor:",ChangePrecision(DEInt`Factor,10);
	print "Precision:",DEInt`Prec;
end intrinsic;


//////////////////////////////////////////////////////////
// ***** Double-exponential numerical integration ***** // 
//////////////////////////////////////////////////////////


function DE_Integrals_Factor( VectorIntegral, Edge, SEC )
// Factors for integrals
	m := SEC`Degree[1];
	n := SEC`Degree[2];
	DFF := SEC`HolomorphicDifferentials;
	ElementaryIntegrals := []; // Needed for Abel-Jacobi
	Integrals := []; // g x (m-1) array of integrals

	// ((b-a)/2)^i, i = 1,..,n
	Fact1 := [ Edge`Data[n-1] ];
	for j in [2..DFF[3][DFF[2]]+1] do
		Append(~Fact1,Edge`Data[n-1] * Fact1[j-1]);
	end for;
	// (2/(b-a))^(nj/m), j = 1,..,m-1
	z := Exp( (n/m) * -Log(Edge`Data[n-1]) ); Fact2 := [z];
	for j in [2..DFF[2]+DFF[1]-1] do
		Append(~Fact2,z*Fact2[j-1]);
	end for;

	// Shift integral back by powers of (p+a)/(p-a)
	ct := 0;
	for j in [1..DFF[2]] do
		PolynomialShiftVector(~VectorIntegral,Edge`Data[n],DFF[3][j],ct);
		ct +:= DFF[3][j];
	end for;

	// Multiply by correct power of zeta, (1-zeta) and ((p-a)/2)^(i+1-jn/m)
	ct := 1; 
	for j in [1..DFF[2]] do
		for ij in [0..DFF[3][j]-1] do
			Pow := -((Edge`up + 1) mod 2)*DFF[4][ct] mod (2*m) + 1;
			VectorIntegral[ct] *:= SEC`Zetas[Pow] * Fact1[ij+1] * Fact2[DFF[4][ct]];
			Append(~ElementaryIntegrals, VectorIntegral[ct]);
			Pow := -2*DFF[4][ct] mod (2*m) + 1;
			VectorIntegral[ct] *:= (1 - SEC`Zetas[Pow]);
			NextIntegrals := [ VectorIntegral[ct] ];
			for l in [1..m-2] do
				Pow := -2*l*DFF[4][ct] mod (2*m) + 1;
				Append(~NextIntegrals,SEC`Zetas[Pow]*VectorIntegral[ct]);
			end for; 
			Append(~Integrals,NextIntegrals);
			ct +:= 1;
		end for;
	end for;
	return Integrals, ElementaryIntegrals;
end function;


function DE_Integrals(Edge,SEC)
	C_0 := Zero(SEC`ComplexField); 
	m := SEC`Degree[1]; 
	nm2 := SEC`Degree[2]-2;
	DFF := SEC`HolomorphicDifferentials;
	DEInt := SEC`IntegrationSchemes[Edge`IntSch];
	VectorIntegral := [ C_0 : j in [1..SEC`Genus] ];

	// Start with zero
	y := DEInt`Weights[1][2]/AC_mthRoot(0,Edge,SEC`Zetas,m,nm2); // 1/y(0)
	wy := y^DFF[1] * DEInt`Weights[1][1];
	ct := 1;
	for j in [1..DFF[2]] do
		if j gt 1 then
			wy *:= y;
		end if;
		wyx := wy;
		VectorIntegral[ct] +:= wyx;
		ct +:= DFF[3][j];
	end for;

	// Evaluate differentials at abisccsas
	for t in [2..DEInt`NPoints] do
		x := DEInt`Abscissas[t];
		mx := -x;
		y1 := DEInt`Weights[t][2]/AC_mthRoot(x,Edge,SEC`Zetas,m,nm2); // 1/y(x)
		y2 := DEInt`Weights[t][2]/AC_mthRoot(mx,Edge,SEC`Zetas,m,nm2); // 1/y(-x)
		wy1 := y1^DFF[1] * DEInt`Weights[t][1];
		wy2 := y2^DFF[1] * DEInt`Weights[t][1];
		ct := 1;
		for j in [1..DFF[2]] do
			if j gt 1 then
				wy1 *:= y1;
				wy2 *:= y2;
			end if;
			wyx1 := wy1;
			wyx2 := wy2;
			VectorIntegral[ct] +:= wyx1 + wyx2;
			ct +:= 1;
			for k in [1..DFF[3][j]-1] do
				wyx1 *:= x;
				wyx2 *:= mx;
				VectorIntegral[ct] +:= wyx1 + wyx2;
				ct +:= 1;
			end for;
		end for;
	end for;
	for j in [1..SEC`Genus] do
		VectorIntegral[j] *:= DEInt`Factor;
	end for;
	return VectorIntegral;
end function;


function DE_Integrals_Edge(EdgeData,SEC)
// Integrate an edge (of the spanning tree) with double-exponential integration }

	// Numerical integration
	VectorIntegral := DE_Integrals(EdgeData,SEC);

	// Multiplication by constants
	PeriodsEdge, ElemIntegralEdge := DE_Integrals_Factor( VectorIntegral,EdgeData,SEC );

	return PeriodsEdge, ElemIntegralEdge;
end function;



function DE_Integrals_Tree(SEC)
// Compute integrals for spanning tree
	Periods := []; ElementaryIntegrals := [];
	for k in [1..SEC`Degree[2]-1] do
		vprint SE,2 : "Integrating edge #:",k;
		P, EI := DE_Integrals_Edge(SEC`SpanningTree`Edges[k],SEC);
		Append(~Periods,P);
		Append(~ElementaryIntegrals,EI);
	end for;
	return Periods, ElementaryIntegrals;
end function;


//////////////////////////////////////////////////////////
// ***** Numerical Integration of Abel-Jacobi map ***** //
//////////////////////////////////////////////////////////


function DE_Integrals_Factor_AJM( VectorIntegral,Edge,SEC)
// EdgeData = [ u_1 , ... , u_{n-1}, (p_x-a)/2, (p_x+a)/(p_x-a), up, <p_x,p_y> 
	C<I> := SEC`ComplexField; 
	g := SEC`Genus; 
	m := SEC`Degree[1];
	n := SEC`Degree[2];
	DFF := SEC`HolomorphicDifferentials;
	// Array of  g integrals
	Integrals := []; 

	// ((p-a)/2)^i, i = 1,..,max_i+1
	Fact1 := [ Edge`Data[n] ];
	for j in [2..DFF[3][DFF[2]]+1] do
		Append(~Fact1,Edge`Data[n] * Fact1[j-1]);
	end for;
	// (2/(p-a))^(nj/m), j = 1,..,m-1
	z := Exp( (n/m) * -Log(Edge`Data[n]) ); Fact2 := [z];
	for j in [2..DFF[2]+DFF[1]-1] do
		Append(~Fact2,z*Fact2[j-1]);
	end for;

	// Shift integral back by powers of (p+a)/(p-a)
	ct := 0;
	for j in [1..DFF[2]] do
		PolynomialShiftVector(~VectorIntegral,Edge`Data[n+1],DFF[3][j],ct);
		ct +:= DFF[3][j];
	end for;

	// Correct sheet?
	P_y_AC := Exp( (n/m) * Log(2*Edge`Data[n]) + (1/m) * Log(SEC`LeadingCoeff))  * SEC`Zetas[Edge`up mod (2*m) + 1]  * AC_mthRoot(1,Edge,SEC`Zetas,m,n-1);

	// Shifting number at P_y
	s := Round((m/(2*Real(Pi(C)))) * ( Arg(Edge`EP[2][2]) - Arg(P_y_AC) ));

	// Multiply by correct power of zeta and ((p-a)/2)^(i+1-jn/m)
	ct := 1;
	for j in [1..DFF[2]] do
		for ij in [0..DFF[3][j]-1] do
			Pow := -(Edge`up + 2*s)*DFF[4][ct] mod (2*m) + 1;
			VectorIntegral[ct] *:= SEC`Zetas[Pow] * Fact1[ij+1] * Fact2[DFF[4][ct]];
			ct +:= 1;
		end for;
	end for;

	return VectorIntegral;
end function;


function DE_Integrals_AJM(Edge,SEC)
	C_0 := Zero(SEC`ComplexField); 
	m := SEC`Degree[1]; 
	nm1 := SEC`Degree[2]-1;
	DFF := SEC`HolomorphicDifferentials;
	VectorIntegral := [ C_0 : j in [1..SEC`Genus] ];
	DEInt := SEC`IntegrationSchemes[Edge`IntSch];

	// Start with zero
	y := DEInt`ExtraFactors[1][1]*DEInt`Weights[1][2]/AC_mthRoot(0,Edge,SEC`Zetas,m,nm1); // 1/y(0)
	wy := y^DFF[1] * DEInt`Weights[1][1];
	ct := 1;
	for j in [1..DFF[2]] do
		if j gt 1 then
			wy *:= y;
		end if;
		wyx := wy;
		VectorIntegral[ct] +:= wyx;
		ct +:= DFF[3][j];
	end for;

	// Evaluate differentials at abisccsas
	for t in [2..DEInt`NPoints] do
		x := DEInt`Abscissas[t];
		mx := -x;
		y1 := DEInt`ExtraFactors[t][1]*DEInt`Weights[t][2]/AC_mthRoot(x,Edge,SEC`Zetas,m,nm1); // 1/y(x)
		y2 := DEInt`ExtraFactors[t][2]*DEInt`Weights[t][2]/AC_mthRoot(mx,Edge,SEC`Zetas,m,nm1); // 1/y(-x)
		wy1 := y1^DFF[1] * DEInt`Weights[t][1];
		wy2 := y2^DFF[1] * DEInt`Weights[t][1];
		ct := 1;
		for j in [1..DFF[2]] do
			if j gt 1 then
				wy1 *:= y1;
				wy2 *:= y2;
			end if;
			wyx1 := wy1;
			wyx2 := wy2;
			VectorIntegral[ct] +:= wyx1 + wyx2;
			ct +:= 1;
			for k in [1..DFF[3][j]-1] do
				wyx1 *:= x;
				wyx2 *:= mx;
				VectorIntegral[ct] +:= wyx1 + wyx2;
				ct +:= 1;
			end for;
		end for;
	end for;
	for j in [1..SEC`Genus] do
		VectorIntegral[j] *:= DEInt`Factor;
	end for;
	return VectorIntegral;
end function;


function DE_Integrals_Edge_AJM(Edge,SEC)
// Integrate an "complex edge" for the Abel-Jacobi map

	// Numerical integration on [-1,1]
	VectorIntegral := DE_Integrals_AJM(Edge,SEC);

	// Multiplication with constants
	Integral_Edge_AJM := DE_Integrals_Factor_AJM(VectorIntegral,Edge,SEC);

	return Integral_Edge_AJM;
end function;
