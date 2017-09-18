/*******************************************************************************

 Copyright (C) 2017 
 Code adapted from Pascal Molin,
 written by Christian Neurohr

*******************************************************************************/


// Import functions;
import "se_spanning_tree.m": SpanningTree, TreeData;
import "se_period_matrix.m": SE_PeriodMatrix;
import "se_abel_jacobi.m": SE_TreeMatrix, SE_ReductionMatrix, SE_RamificationPoints_AJM;
import "se_help_funcs.m": SE_DKPEB, SE_OrdFldComElt;

// Verbose
declare verbose SE,3;

// Superelliptic curves
declare type SECurve;

declare attributes SECurve: 
	DefiningPolynomial, 
	ComplexPolynomial,
	Genus, 
	Degree, 
	BranchPoints, 
	LowPrecBranchPoints, 
	ComplexField, 
	RealField, 
	HolomorphicDifferentials, 
	Prec, 
	SmallPeriodMatrix, 
	BigPeriodMatrix,
	ReductionMatrix, 
	AJM_RamificationPoints, 
	AJM_InftyPoints, 
	AJM_SumOfInftyPoints, 
	TreeMatrix, 
	Error, 
	IntegrationType, 
	InftyPoint,
	ElementaryIntegrals, 
	IntersectionMatrix, 
	Basepoint, 
	Zetas, 
	LeadingCoeff, 
	LeadingCoeff_T, 
	LC2, 
	DifferentialChangeMatrix, 
	SpanningTree;

// Constructor via multivariate polynomial
intrinsic SE_Curve( f::RngMPolElt : Prec := 40, Small := true, AbelJacobi := true, IntegrationType := "Opt", InftyPoints := false ) -> SECurve
{ Creates an superelliptic curve object with defining polynomial f(x,y) = y^m - f(x) }
	K := BaseRing(Parent(f));
	Kx<x> := PolynomialRing(K);
	g := -Evaluate(f,[x,0]);
	m := Degree(f,2);
	return SE_Curve( g, m : Prec := Prec, Small := Small, AbelJacobi := AbelJacobi, IntegrationType := IntegrationType, InftyPoints := InftyPoints );
end intrinsic; 

// Constructor via univariate polynomial and degree
intrinsic SE_Curve( f::RngUPolElt, m::RngIntElt : Prec := 40, Small := true, AbelJacobi := true, IntegrationType := "Opt", InftyPoints := false ) -> SECurve
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
		require Precision(K) ge 40 : "Please enter polynomial with at least 40 digits precision or as exact polynomial.";
		Prec := Min(Prec,Precision(K));
	else
		require K eq Rationals() : "Please enter a polynomial defined over \Q,\R or \C.";
		if Prec lt 40 then
			Prec := 40;
			print "Precision has been increased to 40 decimal digits.";
		end if;
	end if;
	vprint SE,1 : "Precision:",Prec; 

	// Degrees
	n := Degree(f); delta := Gcd(m,n);
	SEC`Degree := [m,n,delta];

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

	// Error
	SEC`Error := 10^-Prec;

	// Estimate lower bound for precision
	fmonic := f/LeadingCoefficient(f);
	if IsExact(K) then
		Fac := Factorization(fmonic);
		Points := []; QPoints := []; MinPrec := 0;
		for k in [1..#Fac] do
			F := Fac[k][1];
			if Degree(F) gt 1 then
				QQ<a> := NumberField(F);
				Cons := Conjugates(a);
				NewPrec := Precision(Universe(Cons));
				if MinPrec eq 0 then
					MinPrec := NewPrec;
				else
					MinPrec := Min(MinPrec,NewPrec);
				end if;
				Points cat:= ChangeUniverse(Cons,ComplexField(NewPrec));
			else
				Append(~QPoints,-Coefficient(F,0));
			end if;
		end for;
		if #Points ne 0 then
			MinPrec := Precision(Universe(Points));
		else
			MinPrec := Prec;
		end if;
		C<I> := ComplexField(MinPrec);
		Points cat:= ChangeUniverse(QPoints,C); 
	else
		CoeffAbs := [Abs(c):c in Coefficients(fmonic) | c ne 0]; 
		MinCH := Floor(Log(10,Min(CoeffAbs)));
		require -MinCH le Prec : "Polynomial coefficients are too small. Please input polynomial to higher precision.";
		MaxCH := Ceiling(Max(0,Ceiling(Log(10,Max(CoeffAbs)))));
		MinPrec := Prec + 2*MaxCH;
		Points := Roots(ChangeRing(fmonic,ComplexField(MinPrec)));
		Points := [ P[1] : P in Points ];
	end if;
	Sort(~Points,SE_OrdFldComElt);
	require #Points eq n : "Defining polynomial has to be squarefree.";

	// Low precision branch points
	SEC`LowPrecBranchPoints := ChangeUniverse(Points,ComplexField(40));
	
	// Increase precision if necessary
	//Prec := Max(Prec,MinimalPrecision);

	// Increase precision for precomputations
	CompPrec := Prec+3;
	SEC`Prec := CompPrec;
	SEC`RealField := RealField(CompPrec);
	vprint SE,1 : "Computational precision:",CompPrec;

	// Compute spanning tree
	SpanningTree(SEC);

	// Maximal M_1
	MaxM1 := Max( [ P[1] : P in SEC`SpanningTree`Params ]);

	// Extra precision
	ExtraPrec := Max(5,Ceiling(Log(10,Binomial(n,Floor(n/4))*MaxM1)));
	vprint SE,1 : "Extra precision:",ExtraPrec; 	

	// Complex field of maximal precision
	MaxPrec := Max(Round(MinPrec/2),CompPrec+ExtraPrec);
	C<I> := ComplexField(MaxPrec);
	if MaxPrec gt MinPrec then
		Points := SE_DKPEB(f,Points,MaxPrec);
	end if;
	SEC`ComplexField := C; 
	vprint SE,1 : "SEC`ComplexField:",SEC`ComplexField;	

	// Embed univariate Polynomial
	f_x := ChangeRing(f,C);
	SEC`ComplexPolynomial := f_x;

	// Leading coefficient
	SEC`LeadingCoeff := LeadingCoefficient(f_x);
	vprint SE,2 : "Leading coefficient:",ChangePrecision(SEC`LeadingCoeff,10); 	

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

	// Transformation, if m = delta (not working properly) // Maybe implemented later
	if false then
	//if n mod m eq 0 then
		CC<I> := ComplexField(2*Precision(SEC`ComplexField));
		SEC`ComplexField := CC;
		D := ScalarMatrix(CC,g,0);
		CCx<x> := PolynomialRing(CC);
		f_x := ChangeRing(f,CC);
		Roots_fx := Roots(f_x);
		Points := [ R[1] : R in Roots_fx ];
		x_n := Points[n];
		f_x := &*[ (x_n - x_k )*x+1 : x_k in Prune(Points)];
		SEC`LeadingCoeff_T := LeadingCoefficient(f_x);
		Roots_f := Roots(f_x);
		Points := [ R[1] : R in Roots_f ];
		require #Points eq n-1 : "f(x) has to be separable.";
		D_LC := [ 1/(SEC`LeadingCoeff^(DFF[4][k]/m) * SEC`LeadingCoeff_T^(DFF[4][k]/m)) : k in [1..g] ];
		print "D_LC:",D_LC;
		SEC`InftyPoint := C!x_n;
		if Abs(x_n) lt SEC`Error then
			D := -DiagonalMatrix(CC,g,D_LC);
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

	// Change due to transformations
	SEC`DifferentialChangeMatrix := D;

	// Branch points
	SEC`BranchPoints := Points;

	// Root of unity powers
	SEC`Zetas := [ Exp(k*Pi(C)*I/m) : k in [0..2*m-1] ];

	// Compute spanning tree
	TreeData(~SEC`SpanningTree,SEC`BranchPoints);

	// Compute big (and small) period matrix
	SE_PeriodMatrix(SEC:Small:=Small,ReductionMatrix:=AbelJacobi);

	// Initiate Abel-Jacobi
	if AbelJacobi then
		// Choose 1st branch point as base point P_0
		SEC`Basepoint := 1;

		// Compute 'map' of the spanning tree
		SE_TreeMatrix(SEC);

		// Compute Abel-Jacobi map between P_0 and other ramification points and sum of infinite points
		SE_RamificationPoints_AJM(SEC);

		SEC`AJM_InftyPoints := [];
		// Compute Abel-Jacobi map between P_0 and P_{\infty}
		if InftyPoints then
			for k in [1..delta] do
				SE_AJM_InftyPoints(k,SEC);
			end for;
			// Test results
			V := (m/delta) * &+[ v : v in SEC`AJM_InftyPoints ];
			assert &and[ Abs(V[j][1]-Round(V[j][1])) lt SEC`Error : j in [1..2*SEC`Genus] ];
		end if;
	end if;
	// Reset to original precision
	SEC`Prec := Prec;
	SEC`RealField := RealField(Prec);
	return SEC;
end intrinsic;


// Constructor via branch points and degree
intrinsic SE_Curve( Points::SeqEnum[FldComElt], m::RngIntElt : LeadingCoeff := 1, Prec := 40, Small := true, AbelJacobi := true, IntegrationType := "Opt", InftyPoints := false) -> SECurve
{ Creates an superelliptic curve object via y^m = LeadingCoeff * \prod[(x-p) : p in Points] }
	Cx<x> := PolynomialRing(Universe(Points));
	f := LeadingCoeff * &*[ (x - p) : p in Points ];
	return SE_Curve(f,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,IntegrationType:=IntegrationType,InftyPoints:=InftyPoints);
end intrinsic; 


// Printing
intrinsic Print( SEC::SECurve : Extended := false )
{ Print data the of complex superelliptic curve SEC }
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
		print " TreeMatrix:", assigned SEC`TreeMatrix;
		print " ReductionMatrix:", assigned SEC`ReductionMatrix;
		print " Basepoint:", SEC`Basepoint;
		print " AJM_RamificationPoints:", assigned SEC`AJM_RamificationPoints;
		print " AJM_SumOfInftyPoints:", assigned SEC`AJM_SumOfInftyPoints;
		print " AJM_InftyPoints:", #SEC`AJM_InftyPoints eq SEC`Degree[3];
	end if;
end intrinsic;

function IsPoint(Point,SEC,Err)
// Checks, whether [x,y] is a point on the curve or not, up to some error
	t := Abs(Point[2]) - Abs(SE_Branches(Point[1],SEC:Principal));
	if Abs(t) gt Err then
		print "SEC:",SEC;
		print "Point:",Point;
		print "t:",t;
		return false;
	else
		return true;
	end if;
end function;


intrinsic SE_RandomCurve( m::RngIntElt, n::RngIntElt : Prec := 40, Ht := 10^2, Monic:=false ) -> RngMPolElt
{ Returns a random superelliptic curve object defined over the rationals }
	Qxy<x,y> := PolynomialRing(Rationals(),2);
	while true do
		if Monic then
			f := x^n;
		else
			f := Random([-1,1])*Random([1..Ht])*x^n;
		end if;
		for j in [0..n-1] do
			c := Random([-1,0,1]);
			f +:= c*Random([1..Ht])*x^j;
		end for;
		g := y^m - f;
		C := Curve(AffineSpace(BaseRing(Parent(g)),2),g);
		if IsAbsolutelyIrreducible(C) and IsSquarefree(UnivariatePolynomial(f)) then
			break;
		end if;
	end while;
	return SE_Curve(g:Prec:=Prec,InftyPoints);
end intrinsic;


// (Simple) Divisor-class on complex superelliptic curves
declare type SEDivisor;

declare attributes SEDivisor: Curve, Support, Degree;

// Constructor
intrinsic SE_Divisor(Points::SeqEnum,SEC::SECurve:Check:=true) -> SEDivisor
{ Creates a divisor on the superelliptic curve SEC }
	// Create object
	SED := New(SEDivisor);	
	SED`Degree := &+[ P[2] : P in Points ];
	SED`Curve := <SEC`Degree[1],SEC`DefiningPolynomial>;

	// Error
	Err := 10^-Floor(SEC`Prec/2);
	//Err := SEC`Error;

	// Complex field
	C<I> :=  SEC`ComplexField;
	C_0 := Zero(C);
	C_1 := One(C);

	// Support and coefficients
	Support := [];
	for P in Points do
		require Type(P[2]) eq RngIntElt and P[2] ne 0 : "Multiplicity (second entry) of",P,"has to be a non-zero integer!";
		if #P[1] eq 1 then
			require Floor(P[1][1]) in [1..SEC`Degree[3]] : "Please enter Infty points as <[k],c_k>, where k in",[1..SEC`Degree[3]]," and c_k is a non-zero integer.";
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
			require P_z ne C_0 : "Please enter Infty points as <[k],c_k>, where k in",[1..SEC`Degree[3]]," and c_k is a non-zero integer.";
			// Embed projective coordinates into affine space over complex field
			P_x := C!P[1][1]/P_z;
			P_y := C!P[1][2]/P_z;
			// Is P a point on the curve?
			if Check then
				require IsPoint([P_x,P_y],SEC,Err) : "Input",P[1],"is not a point on the curve to given precision.";
			end if;
			Append(~Support,<[P_x,P_y],P[2]>);
		else
			require false : "Please input SEDivisors as list of tuples [< P,c_P >,...], where P = [x,y] or P = [x,y,z] and non-zero z for finite points and P = [k] for Infty points";   
		end if;
	end for;
	SED`Support := Support;
	return SED;
end intrinsic;


// Addition
intrinsic '+'(D1::SEDivisor,D2::SEDivisor) -> SEDivisor
{ Add two divisors on SEC }
	require D1`Curve eq D2`Curve : "Divisors have to be defined on the same curve.";
	D1`Support cat:= D2`Support;
	D1`Degree +:= D2`Degree;
	return D1;
end intrinsic;
	

// Printing
intrinsic Print( SED::SEDivisor )
{ Print Riemann surface }
	print " Divisor on superelliptic curve defined as degree",SED`Curve[1],"cover of",SED`Curve[2];
	print " Degree:",SED`Degree;
	print " Support:",SED`Support;
end intrinsic;


intrinsic SE_RandomPoint( SEC::SECurve : Ht := 10^2 ) -> SeqEnum
{ Returns a random point (x,y) on the superelliptic curve SEC }
	m := SEC`Degree[1];
	n := SEC`Degree[2];
	C<I> := ComplexField((Ceiling((m/10))+1)*SEC`Prec+2*n);
	x := C!Random([-Ht..Ht])/Random([1..Ht]) + I*C!Random([-Ht..Ht])/Random([1..Ht]);
	y := SE_Branches(x,SEC:Principal);	
	return [x,y];
end intrinsic;


intrinsic SE_RandomDivisor( NrOfPts::RngIntElt,SEC::SECurve : Ht := 10^2, Zero := false ) -> SEDivisor
{ Creates a random divisor on the superelliptic curve SEC }
	Points := [];
	for k in [1..NrOfPts] do
		Append(~Points,<SE_RandomPoint(SEC:Ht:=Ht),Random([-1,1])*Random([1..Ht])>);
	end for;
	if Zero then
		Deg := &+[ P[2] : P in Points ];
		Points[1][2] -:= Deg;
	end if;
	return SE_Divisor(Points,SEC);
end intrinsic;









