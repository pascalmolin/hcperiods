/******************************************************************************

 Copyright (C) 2017 Christian Neurohr

 ******************************************************************************/

import "se_anal_cont.m": AC_mthRoot;
import "se_help_funcs.m": MakeCCVector, PolynomialShiftVector;
import "se_spanning_tree.m": DE_AJM_Weight;

///////////////////////////////////////////////////////////////
// ***** Parameters for double-exponential integration ***** //
///////////////////////////////////////////////////////////////

C20<i> := ComplexField(20);
RPi20 := Real(Pi(C20));

function Distance_1(P)
	xP := Abs(Re(P));
	i := Parent(P).1;
	if xP gt 1 then
		return Abs(xP-1+i*Im(P));
	else
		return Abs(Im(P));
	end if;
end function;

function Bound_M1(CCV,len,m)
	B1 := 1/(&*[ Distance_1(CCV[k]) : k in [1..len] ])^(1/m);
	if B1 gt 1 then
		return B1^(m-1);
	else
		return B1;
	end if;
end function;

function Distance_thsh( P, r : Lambda := RPi20/2 )
	P := Abs(Re(P)) + i*Abs(Im(P));
	x0 := 0; x1 := Argcosh(RPi20/(2*Lambda*Sin(r))); // s.t. \Lambda Cosh(x)Sin(r)= Pi/2;
	Phi := function(x)
		return Tanh(Lambda*Sinh(x+i*r));
	end function;
	ArgPhiP := function(x)
		return Arg(Cosh(x+i*r))-2*Arg(Cosh(Lambda*Sinh(x+i*r)));
	end function;
	ArgPhi := function(x)
		return Arg(P - Phi(x));
	end function;
	// x := Solve(x := x0, x1, ArgPhi(x)-ArgPhiP(x)-C_Pi/2);
	x := x0;
	Val := ArgPhi(x)-ArgPhiP(x)-RPi20/2;
	n := 10;
	for t in [0..n] do
		x_new := (t/n)*x1 + (1-(t/n))*x0;
		Val_new := ArgPhi(x_new)-ArgPhiP(x_new)-RPi20/2;
		if Abs(Val_new) lt Abs(Val) then
			x := x_new;
			Val := Val_new;
		end if;
	end for;
	return Abs(P-Phi(x));
end function;


function Bound_M2(CCV,len,m,r)
// Compute the bound M2 on the integrand
	B2 := 1/(&*[ Distance_thsh(CCV[k],r) : k in [1..len] ])^(1/m);
	//print "B2:",B2;
	if B2 gt 1 then
		return B2^(m-1);
	else
		return B2;
	end if;
end function;

// CurveData = [n,g,lc_m,C,m]
function DE_Int_Params_Tree( SEC )
// Computes double-exponential integration integration parameters for a spanning tree
	
	C<i> := SEC`ComplexField;
	RPi2 := Real(Pi(C))/2;
	m := SEC`Degree[1];
	n := SEC`Degree[2];
	
	// Compute D
	D := SEC`Prec * Log(C!10);

	// Check r
	r := SEC`SpanningTree`IntPar;
	if r ge RPi2 then
		error Error("Not supposed to happen.");
	end if;
	assert &and[r lt RPi2, r gt 0];

	// Make it precise
	r *:= (19/20);	

	// Parameter Lambda = Pi/2 for DE integration
	Lambda := RPi2;	

	// Compute bounds M1,M2
	M1 := 0; M2 := 0;
	//if #MST gt 0 then
	if true then
		for k in [1..SEC`SpanningTree`Len] do
			Data := SEC`SpanningTree`Data[k];
			LowPrecData := ChangeUniverse(Data,C20);
			M1 := Max(M1,Bound_M1(LowPrecData,n-2,m));
			M2 := Max(M2,Bound_M2(LowPrecData,n-2,m,r));
		end for;
	else
		M1 := 1;
		M2 := 10;
	end if;
	M1 := Ceiling(M1); M2 := Ceiling(M2);
	// New parameters
	Alpha := 1/m;
	//Alpha := (m-1)/m;
	r := C!r; 
	D := C!D;
	X_r := Cos(r) * Sqrt( (Pi(C) / (2*Lambda*Sin(r))) - 1 );
	B_ra := (2/Cos(r)) * ( (X_r/2) * (1/(Cos(Lambda*Sin(r))^(2*Alpha)) + (1/X_r^(2*Alpha))) ) + (1/(2*Alpha*Sinh(X_r)^(2*Alpha)));
	h := Real( 4 * RPi2 * r /  ( D+Log(2*M2 * B_ra + 1)));
	N := Ceiling(Argsinh((D+ Log((2^(2*Alpha+1)*M1)/Alpha )) / ( 2*Alpha*Lambda ))/h);
	
	return DE_Integration(h,N,Lambda,SEC`Degree[1]);
end function;


function DE_Int_Params_AbelJacobi( Edges, SEC )
// Computes double-exponential integration integration parameters for a spanning tree
	
	C<i> := SEC`ComplexField;
	RPi2 := Real(Pi(C))/2;
	m := SEC`Degree[1];
	n := SEC`Degree[2];
	
	// Compute D
	D := SEC`Prec * Log(C!10);

	// Compute r
	r := 4.;
	// Edge = <index,complex number>
	for E in Edges do
		r_E := DE_AJM_Weight(E,SEC`LowPrecBranchPoints,n);
		r := Min(r,r_E);
	end for;
	if r ge RPi2 then
		error Error("Not supposed to happen.");
	end if;
	assert &and[r lt RPi2, r gt 0];

	// Make it precise
	r *:= (19/20);	

	// Parameter Lambda = Pi/2 for DE integration
	Lambda := RPi2;	

	// Compute bounds M1,M2
	M1 := 0; M2 := 0;
	//if #MST gt 0 then
	if true then
		for E in Edges do
			LowPrecData := MakeCCVector(E,SEC`LowPrecBranchPoints);
			//LowPrecData := ChangeUniverse(Data,C20);
			M1 := Max(M1,Bound_M1(LowPrecData,n-1,m));
			M2 := Max(M2,Bound_M2(LowPrecData,n-1,m,r));
		end for;
	else
		M1 := 1;
		M2 := 10;
	end if;
	M1 := Ceiling(M1); M2 := Ceiling(M2);
	// New parameters
	Alpha := 1/m;
	//Alpha := (m-1)/m;
	r := C!r; 
	D := C!D;
	X_r := Cos(r) * Sqrt( (Pi(C) / (2*Lambda*Sin(r))) - 1 );
	B_ra := (2/Cos(r)) * ( (X_r/2) * (1/(Cos(Lambda*Sin(r))^(2*Alpha)) + (1/X_r^(2*Alpha))) ) + (1/(2*Alpha*Sinh(X_r)^(2*Alpha)));
	h := Real( 4 * RPi2 * r /  ( D+Log(2*M2 * B_ra + 1)));
	N := Ceiling(Argsinh((D+ Log((2^(2*Alpha+1)*M1)/Alpha )) / ( 2*Alpha*Lambda ))/h);
	
	return DE_Integration(h,N,Lambda,SEC`Degree[1]);
end function;



/////////////////////////////////////////////////////
// ***** Double-exponential integration type ***** //
/////////////////////////////////////////////////////

// Define DE_Int type
declare type DE_Int;
declare attributes DE_Int: Factor, Steplength, Lambda, Abscissas, Weights, NPoints, Degree;


// Integration parameters
procedure DE_IntegrationPoints( DEInt )
// Computes integration points for double-exponential integration on the interval [-1,1]
// Precision of integration parameters only depens on Precision(h)

	h := DEInt`Steplength;
	R := Parent(h);
	DEInt`Abscissas := [Zero(R)];
	DEInt`Weights := [<One(R),One(R)>];
	
	eh := Exp(h); // e^h
	eh_inv := 1/eh; // e^-h
	ekh := 1; // e^0h
	ekh_inv := 1; // e^-0h

	for k in [1..DEInt`NPoints] do
		ekh *:= eh; // e^kh
		ekh_inv *:= eh_inv; // e^-kh
     		sh := (ekh-ekh_inv)/2; // sinh(kh) = (1/2) * (e^kh - e^-kh)
		ch := (ekh+ekh_inv)/2 ; // cosh(kh) = (e^kh + e^-kh)
      		esh := Exp(DEInt`Lambda*sh); // e^(Lambda*sinh(kh))
      		esh_inv := 1/esh; // e^-(Lambda*sinh(kh))
		chsh := (esh+esh_inv)/2; // cosh(Lambda*sinh(kh)))
		shsh := (esh-esh_inv)/2; // sinh(Lambda*sinh(kh))
		Append(~DEInt`Abscissas,shsh/chsh); // tanh(Lambda*sinh(kh)) =  sinh(Lambda*sinh(kh)) / cosh(Lambda*sinh(kh))
		chsh *:= chsh; // cosh(Lambda*sinh(kh))^2
		Append(~DEInt`Weights,<ch/chsh,Exp( (1/DEInt`Degree) * Log(chsh))>); // <cosh(kh)/cosh(Lambda*sinh(kh))^2,cosh(Lambda*sinh(kh))^(2/m)>
     	end for;
end procedure;


// Constructor
intrinsic DE_Integration( h::FldReElt, N::RngIntElt, Lambda::FldReElt, m::RngIntElt ) -> DE_Int
{ Construct the double exponential integration scheme }

	DEInt := New(DE_Int);
	
	DEInt`Steplength := h;
	DEInt`NPoints := N;
	DEInt`Lambda := Lambda;
	DEInt`Degree := m;

	// Factor for integrals
	DEInt`Factor := DEInt`Lambda * DEInt`Steplength;	

	// Compute weights and abscissas
	DE_IntegrationPoints(DEInt);

	return DEInt;
end intrinsic;


// Printing
intrinsic Print(DEInt::DE_Int)
{ print }
	print "Steplength:",DEInt`Steplength;
	print "Number of abscissas:",DEInt`NPoints;
	print "Lambda:",DEInt`Lambda;
	print "Degree:",DEInt`Degree;
	print "Factor:",DEInt`Factor;
end intrinsic;


function DE_Integrals_Factor( VectorIntegral, EdgeData, SEC)
// EdgeData = [ u_1,...,u_{n-2}, (b-a)/2, (b+a)/(b-a), up ]
	C<i> := SEC`ComplexField; 
	C0 := Zero(C);
	g := SEC`Genus; 
	m := SEC`Degree[1];
	n := SEC`Degree[2];
	DFF := SEC`HolomorphicDifferentials;
	ElementaryIntegrals := []; // Needed for Abel-Jacobi
	Integrals := []; // g x (m-1) array of integrals

	// ((b-a)/2)^i, i = 1,..,n
	Fact1 := [ EdgeData[n-1] ];
	for j in [2..n] do
		Append(~Fact1,EdgeData[n-1] * Fact1[j-1]);
	end for;
	// (2/(b-a))^(nj/m), j = 1,..,m-1
	z := Exp( (n/m) * -Log(EdgeData[n-1]) ); Fact2 := [z];
	for j in [2..m-1] do
		Append(~Fact2,z*Fact2[j-1]);
	end for;

	ct := 0;
	for j in [1..DFF[2]] do
		PolynomialShiftVector(~VectorIntegral,EdgeData[n],DFF[3][j],ct);
		ct +:= DFF[3][j];
	end for;
	assert ct eq g;

	ct := 1; up := Floor(EdgeData[n+1]);
	for j in [1..DFF[2]] do
		for ij in [0..DFF[3][j]-1] do
			Pow := -((up + 1) mod 2)*DFF[4][ct] mod (2*m) + 1;
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
	assert ct eq g+1;
	return Integrals, ElementaryIntegrals;
end function;



function DE_Integrals(EdgeData,SEC,DEInt)
	C<i> := SEC`ComplexField;
	C_0 := Zero(C); g := SEC`Genus; m := SEC`Degree[1]; n := SEC`Degree[2];
	DFF := SEC`HolomorphicDifferentials;
	Integrals := [ C_0 : j in [1..g] ];
	up := Floor(EdgeData[n+1]);

	// Evaluate differentials at abisccsas
	for t in [1..DEInt`NPoints] do
		x := DEInt`Abscissas[t];
		y := DEInt`Weights[t][2]/AC_mthRoot(x,EdgeData,SEC`Zetas,up,m,n-2); // 1/y(x)
		if DFF[1] gt 1 then
			wy := y^DFF[1] * DEInt`Weights[t][1];
		else
			wy := y * DEInt`Weights[t][1];
		end if;
		ct := 1;
		for j in [1..DFF[2]] do
			if j gt 1 then
				wy *:= y;
			end if;
			wyx := wy;
			Integrals[ct] +:= wyx;
			ct +:= 1;
			for k in [1..DFF[3][j]-1] do
				wyx *:= x;
				Integrals[ct] +:= wyx;
				ct +:= 1;
			end for;
		end for;
		assert ct eq g+1;
		if t eq 1 then
			continue;
		end if;
		assert t ne 1;
		x := -x;
		y := DEInt`Weights[t][2]/AC_mthRoot(x,EdgeData,SEC`Zetas,up,m,n-2); // 1/y(x)
		if DFF[1] gt 1 then
			wy := y^DFF[1] * DEInt`Weights[t][1];
		else
			wy := y * DEInt`Weights[t][1];
		end if;
		ct := 1;
		for j in [1..DFF[2]] do
			if j gt 1 then
				wy *:= y;
			end if;
			wyx := wy;
			Integrals[ct] +:= wyx;
			ct +:= 1;
			for k in [1..DFF[3][j]-1] do
				wyx *:= x;
				Integrals[ct] +:= wyx;
				ct +:= 1;
			end for;
		end for;
		assert ct eq g+1;
	end for;
	for j in [1..g] do
		Integrals[j] *:= DEInt`Factor;
	end for;
	return Integrals;
end function;


function DE_Integrals_Edge( EdgeData,SEC,DEInt )
// Integrate an edge (of the spanning tree) with double-exponential integration }

	VectorIntegral := DE_Integrals( EdgeData,SEC,DEInt);

	// Multiply with constants
	PeriodsEdge, ElemIntegralEdge := DE_Integrals_Factor( VectorIntegral,EdgeData,SEC );

	return PeriodsEdge, ElemIntegralEdge;
end function;



function DE_Integrals_Tree( SEC, DEInt )
// Compute integrals for spanning tree
	Periods := []; ElementaryIntegrals := [];
	for k in [1..SEC`SpanningTree`Len] do
		P, EI := DE_Integrals_Edge(SEC`SpanningTree`Data[k],SEC,DEInt);
		Append(~Periods,P);
		Append(~ElementaryIntegrals,EI);
	end for;
	return Periods, ElementaryIntegrals;
end function;





