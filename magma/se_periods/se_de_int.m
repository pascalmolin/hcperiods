/******************************************************************************

 Copyright (C) 2017 
 Code adapted from Pascal Molin,
 written by Christian Neurohr

 ******************************************************************************/

import "se_anal_cont.m": AC_mthRoot;
import "se_help_funcs.m": MakeCCVector, PolynomialShiftVector;
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
	ss := Ceiling(100*tmax);

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
			u := (MDTO/&*[ Dist(z,CCV[k]):k in [1..len] ])^(1/m);
		else
			u := (1/&*[ Dist(z,CCV[k]):k in [1..len] ])^(1/m);
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
declare attributes DE_Int: r, Factor, Steplength, Lambda, Abscissas, Weights, NPoints, Degree, ExtraFactors;


// Integration parameters
procedure DE_IntegrationPoints( DEInt )
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
end procedure;


procedure DE_IntegrationPoints_AJM( DEInt )
// Computes integration points and weights for double-exponential integration on the interval [-1,1]

	h := DEInt`Steplength;
	R := Parent(h);
	mi := 1/DEInt`Degree;
	oh := R!(1/2);
	DEInt`Abscissas := [Zero(R)];
	DEInt`Weights := [<One(R),One(R)>];
	DEInt`ExtraFactors := [<One(R),One(R)>];

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
		thsh := shsh/chsh;
		Append(~DEInt`Abscissas,thsh); // tanh(Lambda*sinh(kh)) =  sinh(Lambda*sinh(kh)) / cosh(Lambda*sinh(kh))
		chsh *:= chsh; // cosh(Lambda*sinh(kh))^2
		// <cosh(kh)/cosh(Lambda*sinh(kh))^2,cosh(Lambda*sinh(kh))^(2/m)>
		Append(~DEInt`Weights,<ch/chsh,chsh^mi>);
		// (1-thsh(kh)^(1/m)
		Append(~DEInt`ExtraFactors,<(1-thsh)^mi,(1+thsh)^mi>); 
     	end for;
end procedure;


// Constructor
intrinsic DE_Integration( Params, SEC : AJM := false ) -> SeqEnum[DE_Int]
{ Construct the double exponential integration scheme }

	// Parameters = < M1,M2,[ r_min,..,r_max ]>
	// r_i < r_{i+1}

	R := RealField(Precision(SEC`ComplexField));
	RPi2 := Pi(R)/2;
	Lambda := RPi2;
	m := SEC`Degree[1];

	// Alpha that provides the biggest 
	Alpha := 1/m;

	// Compute D
	D := SEC`Prec * Log(R!10);
	
	// List of integration schemes
	DE_Integrations := [];

	for P in Params do
		// New scheme
		DEInt := New(DE_Int);
		
		// Get parameters;
		M1 := P[1];
		M2 := P[2];
		r := R!P[3];

		// Compute N,h
		Alpha := 1/m;
		X_r := Cos(r) * Sqrt( ( Pi(R) / (2*Lambda*Sin(r))) - 1 );
		B_ra := (2/Cos(r)) * ( (X_r/2) * (1/(Cos(Lambda*Sin(r))^(2*Alpha)) + (1/X_r^(2*Alpha))) ) + (1/(2*Alpha*Sinh(X_r)^(2*Alpha)));
		h := Real( 4 * RPi2 * r /  ( D+Log(2*M2 * B_ra + 1)));
		N := Ceiling(Argsinh((D+ Log((2^(2*Alpha+1)*M1)/Alpha )) / ( 2*Alpha*Lambda ))/h);

		// Assign attributes
		DEInt`r := ChangePrecision(r,10);
		DEInt`Steplength := h;
		DEInt`NPoints := N;
		DEInt`Lambda := Lambda;
		DEInt`Degree := m;

		// Factor for integrals
		DEInt`Factor := DEInt`Lambda * DEInt`Steplength;	

		// Compute weights and abscissas
		if AJM then
			DE_IntegrationPoints_AJM(DEInt);
		else
			DE_IntegrationPoints(DEInt);

		end if;
		Append(~DE_Integrations,DEInt);
	end for;
	return DE_Integrations;
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
end intrinsic;


//////////////////////////////////////////////////////////
// ***** Double-exponential numerical integration ***** // 
//////////////////////////////////////////////////////////


function DE_Integrals_Factor( VectorIntegral, Edge, SEC )
// EdgeData = [ u_1,...,u_{n-2}, (b-a)/2, (b+a)/(b-a), up ]
	C<I> := SEC`ComplexField; 
	C0 := Zero(C);
	g := SEC`Genus; 
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


function DE_Integrals(Edge,SEC,DEInt)
	C<I> := SEC`ComplexField;
	C_0 := Zero(C); g := SEC`Genus; m := SEC`Degree[1]; n := SEC`Degree[2];
	DFF := SEC`HolomorphicDifferentials;
	VectorIntegral := [ C_0 : j in [1..g] ];

	// Evaluate differentials at abisccsas
	for t in [1..DEInt`NPoints] do
		x := DEInt`Abscissas[t];
		y := DEInt`Weights[t][2]/AC_mthRoot(x,Edge`Data,SEC`Zetas,Edge`up,m,n-2); // 1/y(x)
		wy := y^DFF[1] * DEInt`Weights[t][1];
		ct := 1;
		for j in [1..DFF[2]] do
			if j gt 1 then
				wy *:= y;
			end if;
			wyx := wy;
			VectorIntegral[ct] +:= wyx;
			ct +:= 1;
			for k in [1..DFF[3][j]-1] do
				wyx *:= x;
				VectorIntegral[ct] +:= wyx;
				ct +:= 1;
			end for;
		end for;
		if t eq 1 then
			continue;
		end if;
		x := -x;
		y := DEInt`Weights[t][2]/AC_mthRoot(x,Edge`Data,SEC`Zetas,Edge`up,m,n-2); // 1/y(x)
		wy := y^DFF[1] * DEInt`Weights[t][1];
		ct := 1;
		for j in [1..DFF[2]] do
			if j gt 1 then
				wy *:= y;
			end if;
			wyx := wy;
			VectorIntegral[ct] +:= wyx;
			ct +:= 1;
			for k in [1..DFF[3][j]-1] do
				wyx *:= x;
				VectorIntegral[ct] +:= wyx;
				ct +:= 1;
			end for;
		end for;
	end for;
	for j in [1..g] do
		VectorIntegral[j] *:= DEInt`Factor;
	end for;
	return VectorIntegral;
end function;


function DE_Integrals_Edge( EdgeData,SEC,DEInt )
// Integrate an edge (of the spanning tree) with double-exponential integration }

	// Numerical integration
	VectorIntegral := DE_Integrals( EdgeData,SEC,DEInt);

	// Multiplication by constants
	PeriodsEdge, ElemIntegralEdge := DE_Integrals_Factor( VectorIntegral,EdgeData,SEC );

	return PeriodsEdge, ElemIntegralEdge;
end function;



function DE_Integrals_Tree( SEC, DEInts )
// Compute integrals for spanning tree
	Periods := []; ElementaryIntegrals := [];
	for k in [1..SEC`Degree[2]-1] do
		vprint SE,2 : "Integrating edge #:",k;
		P, EI := DE_Integrals_Edge(SEC`SpanningTree`Edges[k],SEC,DEInts[SEC`SpanningTree`Edges[k]`IntSch]);
		Append(~Periods,P);
		Append(~ElementaryIntegrals,EI);
	end for;
	return Periods, ElementaryIntegrals;
end function;


//////////////////////////////////////////////////////////
// ***** Numerical Integration of Abel-Jacobi map ***** //
//////////////////////////////////////////////////////////


function DE_Integrals_Factor_AJM( VectorIntegral,EdgeData,SEC)
// EdgeData = [ u_1 , ... , u_{n-1}, (p_x-a)/2, (p_x+a)/(p_x-a), up, <p_x,p_y> 
	C<I> := SEC`ComplexField; 
	g := SEC`Genus; 
	m := SEC`Degree[1];
	n := SEC`Degree[2];
	DFF := SEC`HolomorphicDifferentials;
	// Array of  g integrals
	Integrals := []; 

	// ((p-a)/2)^i, i = 1,..,max_i+1
	Fact1 := [ EdgeData[n] ];
	for j in [2..DFF[3][DFF[2]]+1] do
		Append(~Fact1,EdgeData[n] * Fact1[j-1]);
	end for;
	// (2/(p-a))^(nj/m), j = 1,..,m-1
	z := Exp( (n/m) * -Log(EdgeData[n]) ); Fact2 := [z];
	for j in [2..DFF[2]+DFF[1]-1] do
		Append(~Fact2,z*Fact2[j-1]);
	end for;

	// Shift integral back by powers of (p+a)/(p-a)
	ct := 0;
	for j in [1..DFF[2]] do
		PolynomialShiftVector(~VectorIntegral,EdgeData[n+1],DFF[3][j],ct);
		ct +:= DFF[3][j];
	end for;

	// Correct sheet?
	up := Floor(EdgeData[n+2]);
	P_y_AC := Exp( (n/m) * Log(2*EdgeData[n]) + (1/m) * Log(SEC`LeadingCoeff))  * SEC`Zetas[up mod (2*m) + 1]  * AC_mthRoot(1,EdgeData,SEC`Zetas,up,m,n-1);

	// Shifting number at P_y
	s := Round((m/(2*Real(Pi(C)))) * ( Arg(EdgeData[n+3]) - Arg(P_y_AC) ));

	// Multiply by correct power of zeta and ((p-a)/2)^(i+1-jn/m)
	ct := 1;
	for j in [1..DFF[2]] do
		for ij in [0..DFF[3][j]-1] do
			Pow := -(up + 2*s)*DFF[4][ct] mod (2*m) + 1;
			VectorIntegral[ct] *:= SEC`Zetas[Pow] * Fact1[ij+1] * Fact2[DFF[4][ct]];
			ct +:= 1;
		end for;
	end for;

	return VectorIntegral;
end function;


function DE_Integrals_AJM( EdgeData, SEC, DEInt )
	C<I> := SEC`ComplexField;
	C_0 := Zero(C); 
	g := SEC`Genus;
	m := SEC`Degree[1]; 
	n := SEC`Degree[2];
	DFF := SEC`HolomorphicDifferentials;
	VectorIntegral := [ C_0 : j in [1..g] ];
	up := Floor(EdgeData[n+2]);

	// Evaluate differentials at abisccsas
	for t in [1..DEInt`NPoints] do
		x := DEInt`Abscissas[t];
		y := DEInt`ExtraFactors[t][1]*DEInt`Weights[t][2]/AC_mthRoot(x,EdgeData,SEC`Zetas,up,m,n-1); // 1/y(x)
		wy := y^DFF[1] * DEInt`Weights[t][1];
		ct := 1;
		for j in [1..DFF[2]] do
			if j gt 1 then
				wy *:= y;
			end if;
			wyx := wy;
			VectorIntegral[ct] +:= wyx;
			ct +:= 1;
			for k in [1..DFF[3][j]-1] do
				wyx *:= x;
				VectorIntegral[ct] +:= wyx;
				ct +:= 1;
			end for;
		end for;
		if t eq 1 then
			continue;
		end if;
		x := -x;
		y := DEInt`ExtraFactors[t][2]*DEInt`Weights[t][2]/AC_mthRoot(x,EdgeData,SEC`Zetas,up,m,n-1); // 1/y(x)
		wy := y^DFF[1] * DEInt`Weights[t][1];
		ct := 1;
		for j in [1..DFF[2]] do
			if j gt 1 then
				wy *:= y;
			end if;
			wyx := wy;
			VectorIntegral[ct] +:= wyx;
			ct +:= 1;
			for k in [1..DFF[3][j]-1] do
				wyx *:= x;
				VectorIntegral[ct] +:= wyx;
				ct +:= 1;
			end for;
		end for;
	end for;
	for j in [1..g] do
		VectorIntegral[j] *:= DEInt`Factor;
	end for;
	return VectorIntegral;
end function;


function DE_Integrals_Edge_AJM( ComplexEdge, SEC, DEInt)
// Integrate an "complex edge" for the Abel-Jacobi map
// ComplexEdge = <P=[x,y],v_P,k> , such that P_k is the closest ramification points

	// Compute data for integration
	EdgeData, up := MakeCCVector(<ComplexEdge[3],ComplexEdge[1][1]>,SEC`BranchPoints);
	Append(~EdgeData,up);
	Append(~EdgeData,ComplexEdge[1][2]);

	// Numerical integration on [-1,1]
	VectorIntegral := DE_Integrals_AJM( EdgeData, SEC, DEInt );

	// Multiplication with constants
	Integral_Edge_AJM := DE_Integrals_Factor_AJM( VectorIntegral,EdgeData,SEC);

	return Integral_Edge_AJM;
end function;
