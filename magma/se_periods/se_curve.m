/******************************************************************************

 Copyright (C) 2017 Christian Neurohr

 ******************************************************************************/

// Import global settings
import "superelliptic.m": RS_SEIntegrate, RS_SEIntegrate3, RS_NRootAC, RS_SEMST, RS_SEMST2, RS_MakeCCVectors, RS_SETau3, RS_SEInfTau;
import "se_spanning_tree.m": SpanningTree, TreeData;
import "se_period_matrix.m": SE_PeriodMatrix;


// Verbose
declare verbose SE,3;

C_20<i> := ComplexField(20);
C_Pi := Real(Pi(C_20));

// Superelliptic curves
declare type SECurve;

declare attributes SECurve: 
	DefiningPolynomial, Genus, Degree, BranchPoints, LowPrecBranchPoints, ComplexField, RealField, HolomorphicDifferentials, Prec, SmallPeriodMatrix, BigPeriodMatrix, ReductionMatrix, 
	AJM_RamPoints, AJM_InftyPoints, AbelJacobi, TreeMatrix, Error, IntegrationType, InfinitePoint,
ElementaryIntegrals, spsm_Matrix, IntersectionMatrix, Basepoint, Zetas, LeadingCoeff, LeadingCoeff_T, LC2, ComplexPolynomial,
DifferentialChangeMatrix, SpanningTree;

// Constructor via multivariate polynomial
intrinsic SE_Curve( f::RngMPolElt : Prec := 20, Small := true, AbelJacobi := false ) -> SECurve
{ Creates an superelliptic curve object with defining polynomial f(x,y) = y^m - f(x) }
	K := BaseRing(Parent(f));
	Kx<x> := PolynomialRing(K);
	g := Evaluate(f,[x,0]);
	m := Degree(f,2);
	return SE_Curve( g, m : Prec := Prec, Small := Small, AbelJacobi := AbelJacobi );
end intrinsic; 

// Constructor via univariate polynomial and degree
intrinsic SE_Curve( f::RngUPolElt, m::RngIntElt : Prec := 20, Small := true, AbelJacobi := false ) -> SECurve
{ Creates an superelliptic curve object with defining equation y^m = f(x) }

	// Create object
	SEC := New(SECurve);	

	// Defining polynomial
	SEC`DefiningPolynomial := f;

	// Precision
	K := BaseRing(Parent(f));
	if not IsExact(K) then
		Prec := Max(Prec,Precision(K));
	end if;

	// Degrees
	n := Degree(f); delta := Gcd(m,n);
	
	// Requirements
	require &and[m ge 2, n ge 3] : "Degrees not supported."; 

	// Extra precision
	ExtraPrec := Max(2*Ceiling(Log(10,Binomial(n,Floor(n/2)))),5);
	vprint SE,1 : "Extra precision:",ExtraPrec; 	

	// Fields
	C20<i> := ComplexField(20);
	C<i> := ComplexField(Prec+ExtraPrec);
	Cx<x> := PolynomialRing(C);
	SEC`ComplexField := C; 
	SEC`RealField := RealField(Prec+5); 
	SEC`Prec := Prec;

	// Error
	SEC`Error := 10^-SEC`Prec;

	// Embed univariate Polynomial
	f_x := ChangeRing(f,C);
	SEC`ComplexPolynomial := f_x;

	// Roots of f_x
	Roots_fx := Roots(f_x);
	Points := [ R[1] : R in Roots_fx ];
	require #Points eq n : "Defining polynomial has to be squarefree.";

	// Leading coefficient
	SEC`LeadingCoeff := LeadingCoefficient(f_x);
	vprint SE,2 : "Leading coefficient:",SEC`LeadingCoeff; 	

	// Holomorphic differentials in the form
	// DFF = [ Min_j, #js, [ Max_ij, j = Min_j to Max_j ], [ jPows ] ];
	jm := Max(1,Ceiling((m+delta)/n));
	DF := []; jPows := [];
	j := jm;
	while j lt m do
		dum := Floor((j*n-delta)/m);
		if dum le n-1 then
			Append(~DF,dum);
			jPows cat:= [ j : k in [1..dum] ];
			j +:= 1;
		else
			break;
		end if;
	end while;
	assert j-jm eq #DF;
	DFF := <jm,#DF,DF,jPows>;
	SEC`HolomorphicDifferentials := DFF;

	// Genus
	g := Round((1/2) * ((m-1)*(n-1) - delta + 1));
	SEC`Genus := g;

	// Transform, if m = delta
	if false then
	//if n mod m eq 0 then
		"TRAAAAAAAAAAAAAAAAAAAANSFORMATION!";
		D := ScalarMatrix(C,g,0); print "D:",D;
		x_n := Points[n];
		f_x := &*[ (x_n - x_k )*x+1 : x_k in Prune(Points)];
		print "f_x:",f_x;
		SEC`LeadingCoeff_T := LeadingCoefficient(f_x);
		print "LC-org:",SEC`LeadingCoeff;
		print "LC2:",SEC`LeadingCoeff_T;
		Roots_f := Roots(f_x);
		Points := [ R[1] : R in Roots_f ];
		print "Points:",Points;
		require #Points eq n-1 : "f(x,0) has to be separable.";
		//LLCS :=  -(1/m) * (Log(SEC`LeadingCoeff) + Log(SEC`LeadingCoeff_T));
		//print "LLCS:",LLCS;
		//D_LC := [ Exp( DFF[j][2] * LLCS ) : j in [1..SEC`Genus] ];
		D_LC := [ 1/(SEC`LeadingCoeff^(DFF[4][k]/m) * SEC`LeadingCoeff_T^(DFF[4][k]/m)) : k in [1..g] ];
		print "D_LC:",D_LC;
		SEC`InfinitePoint := C!x_n;
		print "x_n:",x_n;
		print "Pxn:",Parent(x_n);
		print "C:",C;
		if Abs(x_n) lt SEC`Error then
			"WEIRD CASE";
			D := -DiagonalMatrix(C,g,D_LC);
		else
			M := n/m; 	
			//for j in [1..SEC`Genus] do
			ct := 1;
			print "DFF:",DFF;
			for j in [0..DFF[2]-1] do
 				for ij in [0..DFF[3][j+1]-1] do
					print "CT:",ct;
					nj := DFF[1]+j;
					e := Round(nj*M - ij - 2);
					print "e:",e;
					assert e ge 0;
					for k in [0..e] do
						print "Binomial(e,k) * x_n^(e-k):",-D_LC[nj] * Binomial(e,k) * (-x_n)^(e-k);
						D[k+1][ct] -:= D_LC[nj] * Binomial(e,k) * (-x_n)^(e-k);
					end for;
					ct +:= 1;
				end for;
			end for;
			assert ct eq g+1;
		end if;
		n -:= 1;
		delta := 1;
	else
		LLC1 := -(1/m)*Log(SEC`LeadingCoeff);
		D := DiagonalMatrix(C,SEC`Genus,[ Exp( DFF[4][k] * LLC1 ) : k in [1..SEC`Genus] ]);
	end if;

	// Degrees
	SEC`Degree := [m,n,delta];

	// Change due to transformations
	SEC`DifferentialChangeMatrix := D;

	// Branch points
	SEC`BranchPoints := Points;
	// Low precision branch Points
	SEC`LowPrecBranchPoints := ChangeUniverse(SEC`BranchPoints,C20);

	// Integration method
	if m eq 2 then
		// Gauss-Chebychev
		SEC`IntegrationType := "GC";
	else
		// Double-exponential
		SEC`IntegrationType := "DE";
	end if;
	// Root of unity powers
	SEC`Zetas := [ Exp(k*Pi(C)*i/m) : k in [0..2*m-1] ];

	// Compute spanning tree
	SEC`SpanningTree := SpanningTree(SEC);
	TreeData(~SEC`SpanningTree,SEC`BranchPoints);

	// Choose 1st branch point as base point
	SEC`Basepoint := 1;

	// Compute period matrix
	SE_PeriodMatrix(SEC:Small:=Small);

	// Set up Abel-Jacobi map
	if AbelJacobi then
		SE_InitiateAbelJacobi(SEC:Small:=Small);
	end if;

	// Return curve
	return SEC;
end intrinsic;


// Constructor via branch points and degree
intrinsic SE_Curve( Points::SeqEnum[FldComElt], m::RngIntElt : LeadingCoeff := 1, Prec := 20, Small := true, AbelJacobi := false ) -> SECurve
{ Creates an superelliptic curve object via y^m = LC * prod(p in Points) (x-p)  }
	Cx<x> := PolynomialRing(Universe(Points));
	f := LeadingCoeff * &*[ (x - p) : p in Points ];
	return SE_Curve(f,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi);
end intrinsic; 


// Printing
intrinsic Print( SEC::SECurve )
{ Print Riemann surface }
	print "Superelliptic curve of genus", SEC`Genus ,"defined by: 0 =",SEC`DefiningPolynomial,"and prescribed precision",SEC`Prec;
	print "";
	print "Computed data:";
	print " ComplexField:", SEC`ComplexField;
	print " RealField:", SEC`RealField;
	print " HolomorphicDifferentials:", assigned SEC`HolomorphicDifferentials;
	print " BranchPoints:", assigned SEC`BranchPoints;
	//print " Roots of unity:", assigned SEC`Zetas;
	print " SpanningTree:", assigned SEC`SpanningTree;
	print " IntegrationType:", SEC`IntegrationType;
	//print " spsm_Matrix:", assigned SEC`spsm_Matrix;
	print " IntersectionMatrix:", assigned SEC`IntersectionMatrix;
	print " ElementaryIntegrals:", assigned SEC`ElementaryIntegrals;
	print " BigPeriodMatrix:",assigned SEC`BigPeriodMatrix;
	print " SmallPeriodMatrix:",assigned SEC`SmallPeriodMatrix;
	print " AbelJacobi:", assigned SEC`AbelJacobi;
	print " AJM_RamPoints:", assigned SEC`AJM_RamPoints;
	print " AJM_InftyPoints:", assigned SEC`AJM_InftyPoints;
	print " Basepoint:", assigned SEC`Basepoint;
	print " TreeMatrix:", assigned SEC`TreeMatrix;
	print " ReductionMatrix:", assigned SEC`ReductionMatrix;
end intrinsic;


// Divisors
declare type SEDivisor;

declare attributes SEDivisor: DefiningPolynomial, Support, Coefficients, Degree;

// Constructor
intrinsic SE_Divisor(Points::Tup,Coeffs::SeqEnum[RngIntElt],SEC::SECurve : Check:=true) -> SEDivisor
{ Construct the divisor \sum v_P P }
	require &and[#Points eq #Coeffs,0 notin Coeffs] : "Every point in the support needs a non-zero coefficient.";
	SED := New(SEDivisor);		
	SED`Degree := &+Coeffs;
	SED`Coefficients := Coeffs;	
	SED`DefiningPolynomial := SEC`DefiningPolynomial;

	// Complex field
	C<i> :=  SEC`ComplexField; C_1 := One(C); C_0 := Zero(C);

	// Support
	Supp := < >;
	for P in Points do
		if #P eq 1 then
			x_P := C!P[1]; y_P := SE_PrincipalBranch(SEC,x_P)[1]; 
			Append(~Supp,[x_P,y_P]);
		elif #P eq 2 then
			x_P := C!P[1]; y_P := C!P[2];
			if Check then
				t := Evaluate(SEC`ComplexPolynomial,x_P)- y_P^SEC`Degree[1];
				if Abs(t) gt 10^-10 then
					print "t:",t;
					print "x_P:",x_P; 
					print "y_P:",y_P;
					error Error("Not a point on the curve.");
				end if;
				
			end if;
			Append(~Supp,[C!P[1],C!P[2]]);
		elif #P eq 3 then
			z_P := C!P[3];
			if z_P eq C_0 then
				require &and[ P[1] eq 0,P[2] in [1..SEC`Degree[3]] ] :
					"Please enter points above infinity as [0,j,0], indexed by j in",[1..SEC`Degree[3]];
				Append(~Supp,[0,P[2],0]);
			else
				x_P := P[1]/z_P;
				y_P := P[2]/z_P;
				if Check then
					tt := Evaluate(SEC`ComplexPolynomial,x_P)- y_P^SEC`Degree[1];
					if Abs(t) gt 10^-10 then
						print "t:",t;
						print "x_P:",x_P; 
						print "y_P:",y_P;
						error Error("Not a point on the curve.");
					end if;
				end if;
				Append(~Supp,[C!P[1],C!P[2]]);
			end if;
		else
			require false : "Correct input: Points = < Type1,Type2,Type3,... >, where Type1 = [x_P] (just x-coordinate), Type2 = [x_P,y_P] (point on the curve), Type3 = [0,j,0] (point above infinity)"; 
		end if;
	end for;
	SED`Support := Supp;
	return SED;
end intrinsic;

// Printing
intrinsic Print( SED::SEDivisor )
{ Print Riemann surface }
	print " Divisor on curve given by 0 =",SED`DefiningPolynomial;
	print " Degree:",SED`Degree;
	print " Support:",SED`Support;
	print " Coefficients:",SED`Coefficients; 
end intrinsic;

intrinsic SE_ZeroDivisor(Points::Tup,Coeffs::SeqEnum[RngIntElt],SEC::SECurve) -> SEDivisor
{ Creates a zero divisor on SEC by adding -Deg * P_0 }
	Deg := &+Coeffs;
	Append(~Points,[SEC`BranchPoints[SEC`Basepoint],0]);
	Append(~Coeffs,-Deg);
	return SE_Divisor(Points,Coeffs,SEC);
end intrinsic;	









