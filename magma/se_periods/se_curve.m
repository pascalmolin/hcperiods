/*******************************************************************************

 Copyright (C) 2017 
 Code adapted from Pascal Molin,
 written by Christian Neurohr

*******************************************************************************/


// Import functions;
import "se_spanning_tree.m": SpanningTree, TreeData;
import "se_period_matrix.m": SE_PeriodMatrix;
import "se_abel_jacobi.m": SE_TreeMatrix, SE_ReductionMatrix, SE_RamificationPoints_AJM;

// Verbose
declare verbose SE,3;

C_20<i> := ComplexField(20);
C_Pi := Real(Pi(C_20));

// Superelliptic curves
declare type SECurve;

declare attributes SECurve: 
	DefiningPolynomial, Genus, Degree, BranchPoints, LowPrecBranchPoints, ComplexField, RealField, HolomorphicDifferentials, Prec, SmallPeriodMatrix, BigPeriodMatrix, ReductionMatrix, 
	AJM_RamificationPoints, AJM_InftyPoints, AJM_SumOfInftyPoints, AbelJacobiMap, TreeMatrix, Error, IntegrationType, InfinitePoint,
ElementaryIntegrals, IntersectionMatrix, Basepoint, Zetas, LeadingCoeff, LeadingCoeff_T, LC2, ComplexPolynomial,
DifferentialChangeMatrix, SpanningTree;

// Constructor via multivariate polynomial
intrinsic SE_Curve( f::RngMPolElt : Prec := 20, Small := false, AbelJacobi := true, IntegrationType := "Opt", InfinitePoints := false ) -> SECurve
{ Creates an superelliptic curve object with defining polynomial f(x,y) = y^m - f(x) }
	K := BaseRing(Parent(f));
	Kx<x> := PolynomialRing(K);
	g := -Evaluate(f,[x,0]);
	m := Degree(f,2);
	return SE_Curve( g, m : Prec := Prec, Small := Small, AbelJacobi := AbelJacobi, IntegrationType := IntegrationType, InfinitePoints := InfinitePoints );
end intrinsic; 

// Constructor via univariate polynomial and degree
intrinsic SE_Curve( f::RngUPolElt, m::RngIntElt : Prec := 20, Small := false, AbelJacobi := true, IntegrationType := "Opt", InfinitePoints := false ) -> SECurve
{ Creates an superelliptic curve object with defining equation y^m = f(x) }

	// Create object
	SEC := New(SECurve);	

	// Defining polynomial
	SEC`DefiningPolynomial := f;
	vprint SE,1 : "DefiningPolynomial:",f; 	
	vprint SE,1 : "Degree m:",m; 	

	// Precision
	K := BaseRing(Parent(f));
	if not IsExact(K) then
		Prec := Max(Prec,Precision(K));
	end if;

	// Degrees
	n := Degree(f); delta := Gcd(m,n);

	// Requirements
	require &and[m ge 2, n ge 3] : "Degrees not supported."; 

	// Genus
	g := Round((1/2) * ((m-1)*(n-1) - delta + 1));
	SEC`Genus := g;

	// Integration method
	require IntegrationType in ["Opt","DE","GC"] : "Invalid integration type.";
	if m gt 2 or IntegrationType eq "DE" then
		// Double-exponential
		SEC`IntegrationType := "DE";
	else
		// Gauss-Chebychev
		SEC`IntegrationType := "GC";
	end if;


	// Extra precision
	ExtraPrec := 2*Ceiling(Log(10,Binomial(n,Floor(n/2))))+5;
	vprint SE,1 : "Extra precision:",ExtraPrec; 	

	// Fields
	C20<i> := ComplexField(20);

	// Error
	SEC`Error := 10^-Prec;
	
	if AbelJacobi then
		// Increase precision for Abel-Jacobi map precomputations
		Prec +:= 5;
	end if;

	C<i> := ComplexField(Prec+ExtraPrec);
	Cx<x> := PolynomialRing(C);
	SEC`ComplexField := C; 
	SEC`Prec := Prec;
	SEC`RealField := RealField(Prec+5);
	vprint SE,1 : "Precision:",Prec; 	
	

	// Embed univariate Polynomial
	f_x := ChangeRing(f,C);
	SEC`ComplexPolynomial := f_x;

	// Roots of f_x
	Roots_fx := Roots(f_x);
	Points := [ R[1] : R in Roots_fx ];
	require #Points eq n : "Defining polynomial has to be squarefree.";

	// Leading coefficient
	SEC`LeadingCoeff := LeadingCoefficient(f_x);
	vprint SE,2 : "Leading coefficient:",C20!SEC`LeadingCoeff; 	

	// Holomorphic differentials in the form
	// DFF = [ min_j, #j's, [ max_ij : j in [min_j..max_j] ], [ jPows ] ];
	jm := Max(1,Ceiling((m+delta)/n));
	DF := []; jPows := [];
	j := jm;
	while j lt m do
		k := Floor((j*n-delta)/m);
		if k le n-1 then
			Append(~DF,k);
			jPows cat:= [ j : l in [1..k] ];
			j +:= 1;
		else
			break;
		end if;
	end while;
	assert j-jm eq #DF;
	DFF := <jm,#DF,DF,jPows>;
	SEC`HolomorphicDifferentials := DFF;

	// Transformation, if m = delta (?)
	if false then
	//if n mod m eq 0 then
		D := ScalarMatrix(C,g,0);
		x_n := Points[n];
		f_x := &*[ (x_n - x_k )*x+1 : x_k in Prune(Points)];
		print "f_x:",f_x;
		SEC`LeadingCoeff_T := LeadingCoefficient(f_x);
		print "LC-org:",SEC`LeadingCoeff;
		print "LC2:",SEC`LeadingCoeff_T;
		Roots_f := Roots(f_x);
		Points := [ R[1] : R in Roots_f ];
		print "Points:",Points;
		require #Points eq n-1 : "f(x) has to be separable.";
		D_LC := [ 1/(SEC`LeadingCoeff^(DFF[4][k]/m) * SEC`LeadingCoeff_T^(DFF[4][k]/m)) : k in [1..g] ];
		print "D_LC:",D_LC;
		SEC`InfinitePoint := C!x_n;
		if Abs(x_n) lt SEC`Error then
			D := -DiagonalMatrix(C,g,D_LC);
		else
			M := n/m; 	
			ct := 1;
			for j in [0..DFF[2]-1] do
 				for ij in [0..DFF[3][j+1]-1] do
					nj := DFF[1]+j;
					e := Round(nj*M - ij - 2);
					assert e ge 0;
					for k in [0..e] do
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

	// Root of unity powers
	SEC`Zetas := [ Exp(k*Pi(C)*i/m) : k in [0..2*m-1] ];

	// Compute spanning tree
	SpanningTree(SEC);
	TreeData(~SEC`SpanningTree,SEC`BranchPoints);

	// Compute big (and small) period matrix
	SE_PeriodMatrix(SEC:Small:=Small);

	// Initiate Abel-Jacobi
	if AbelJacobi then
		// Choose 1st branch point as base point P_0
		SEC`Basepoint := 1;

		// Compute 'map' of the spanning tree
		SE_TreeMatrix(SEC);

		// Compute reduction matrix
		SE_ReductionMatrix(SEC);

		// Compute Abel-Jacobi map between P_0 and other ramification points and sum of infinite points
		SE_RamificationPoints_AJM(SEC);

		SEC`AJM_InftyPoints := [];
		// Compute Abel-Jacobi map between P_0 and P_{\infty}
		if InfinitePoints then
			for k in [1..delta] do
				SE_InfinitePoints_AJM(k,SEC);
			end for;
			// Test results
			V := (m/delta) * &+[ v : v in SEC`AJM_InftyPoints ];
			assert &and[ Abs(V[j][1]-Round(V[j][1])) lt SEC`Error : j in [1..2*SEC`Genus] ];
		end if;

		// Reset to original precision
		SEC`RealField := RealField(SEC`Prec);
		SEC`Prec -:= 5;
	end if;
	return SEC;
end intrinsic;


// Constructor via branch points and degree
intrinsic SE_Curve( Points::SeqEnum[FldComElt], m::RngIntElt : LeadingCoeff := 1, Prec := 20, Small := true, AbelJacobi := true, IntegrationType := "Opt", InfinitePoints := false) -> SECurve
{ Creates an superelliptic curve object via y^m = LeadingCoeff * &*[ (x-p) : p in Points] }
	Cx<x> := PolynomialRing(Universe(Points));
	f := LeadingCoeff * &*[ (x - p) : p in Points ];
	return SE_Curve(f,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,IntegrationType := IntegrationType, InfinitePoints := InfinitePoints);
end intrinsic; 


// Printing
intrinsic Print( SEC::SECurve : Extended := false )
{ Print Riemann surface }
	print "Complex superelliptic curve of genus",SEC`Genus ,"defined as degree",SEC`Degree[1],"cover of",SEC`DefiningPolynomial,"with prescribed precision",SEC`Prec;
	if Extended then
		print "";
		print "Computed data:";
		print " ComplexField:", SEC`ComplexField;
		print " RealField:", SEC`RealField;
		print " HolomorphicDifferentials:", assigned SEC`HolomorphicDifferentials;
		print " BranchPoints:", assigned SEC`BranchPoints;
		print " SpanningTree:", assigned SEC`SpanningTree;
		print " IntegrationType:", SEC`IntegrationType;
		print " IntersectionMatrix:", assigned SEC`IntersectionMatrix;
		print " ElementaryIntegrals:", assigned SEC`ElementaryIntegrals;
		print " BigPeriodMatrix:",assigned SEC`BigPeriodMatrix;
		print " SmallPeriodMatrix:",assigned SEC`SmallPeriodMatrix;
		print " Basepoint:", assigned SEC`Basepoint;
		print " TreeMatrix:", assigned SEC`TreeMatrix;
		print " ReductionMatrix:", assigned SEC`ReductionMatrix;
		print " AJM_RamificationPoints:", assigned SEC`AJM_RamificationPoints;
		print " AJM_SumOfInftyPoints:", assigned SEC`AJM_SumOfInftyPoints;
		print " AJM_InftyPoints:", #SEC`AJM_InftyPoints eq SEC`Degree[1];
	end if;
end intrinsic;




function IsPoint(Point,SEC,Err)
// Checks, whether [x,y] is a point on the curve or not, up to some error
	t := Evaluate(SEC`ComplexPolynomial,Point[1])- Point[2]^SEC`Degree[1];
	if Abs(t) gt Err then
		print "SEC:",SEC;
		print "Tup:",Tup;
		print "t:",t;
		return false;
	else
		return true;
	end if;
end function;

// Divisors
declare type SEDivisor;

declare attributes SEDivisor: Curve, Support, Degree;


// Constructor
intrinsic SE_Divisor(Points::SeqEnum[Tup],SEC::SECurve:Check:=true) -> SEDivisor
{ Construct the divisor \sum v_P P }
	SED := New(SEDivisor);	
	SED`Degree := &+[ P[2] : P in Points ];
	SED`Curve := <SEC`Degree[1],SEC`DefiningPolynomial>;
	
	// Error
	Err := 10^-10;

	// Complex field
	C<i> :=  SEC`ComplexField; 
	C_1 := One(C); C_0 := Zero(C);

	// Support
	Support := [];
	for P in Points do
		require Type(P[2]) eq RngIntElt and P[2] ne 0 : "Multiplicity (second entry) of",P,"has to be non-zero integer!";
		if #P[1] eq 1 then
			require Floor(P[1][1]) in [1..SEC`Degree[3]] : "Please enter infinite points as <[k],c_k>, where k in",[1..SEC`Degree[3]]," and c_k is an integer.";
			Append(~Support,P);
		elif #P[1] eq 2 then
			// Embed affine coordinates into complex field
			P_x := C!P[1][1]; 
			P_y := C!P[1][2];
			// Is P a point on the curve?
			if Check then
				require IsPoint([P_x,P_y],SEC,Err) : "Input",P[1],"is not a point on the curve to given precision.";
			end if;
			Append(~Support,<[P_x,P_y],P[2]>);
		elif #P[1] eq 3 then
			P_z := C!P[1][3];
			require P_z ne C_0 : "Please enter infinite points as <[k],c_k>, where k in",[1..SEC`Degree[3]]," and c_k is an integer.";
			// Embed projective coordinates into affine space over complex field
			P_x := C!P[1][1]/P_z;
			P_y := C!P[1][2]/P_z;
			// Is P a point on the curve?
			if Check then
				require IsPoint([P_x,P_y],SEC,Err) : "Input",P[1],"is not a point on the curve to given precision.";
			end if;
			Append(~Support,<[P_x,P_y],P[2]>);
		else
			require false : "Please input SEDivisors as list of tuples [< P,c_P >,...], where P = [x,y] or P = [x,y,z] and non-zero z for finite points and P = [k] for infinite points";   
		end if;
	end for;
	SED`Support := Support;
	return SED;
end intrinsic;


// Printing
intrinsic Print( SED::SEDivisor )
{ Print Riemann surface }
	print " Divisor on superelliptic curve defined as degree",SED`Curve[1],"cover of",SED`Curve[2];
	print " Degree:",SED`Degree;
	print " Support:",SED`Support;
end intrinsic;


intrinsic SE_RandomDivisor(NrOfPts::RngIntElt,SEC::SECurve : Ht := 1000 ) -> SEDivisor
{ Creates a random divisor on the superelliptic curve SEC }
	C<i> := SEC`ComplexField;
	Points := [];
	for k in [1..NrOfPts] do
		x := C!Random([-Ht..Ht])/Random([1..Ht]) + i*C!Random([-Ht..Ht])/Random([1..Ht]);
		y := SE_PrincipalBranch(x,SEC);	
		Append(~Points,<[x,y],Random([-1,1])*Random([1..Ht])>);
	end for;
	return SE_Divisor(Points,SEC);
end intrinsic;









