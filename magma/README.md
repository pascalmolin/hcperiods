///// ***** Quick guide on how to use the MAGMA-package se_periods ***** /////

	Contact: neurohrchristian@googlemail.com

	1.) Copy folder se_periods and spec-file into some directory PATH
	2.) Load MAGMA in terminal
	3.) Load package by typing:	AttachSpec("PATH/spec");
	4.) Use commands in MAGMA, e.g.

		// Compute the period matrices of y^m = f(x) to precision Prec digits (defined over rationals)

		Qx<x> := PolynomialRing(Rationals());
		f := x^5 + x^3 + x + 1;
		m := 3;
		Prec := 200;
		M_1 := SE_BigPeriodMatrix(f,m:Prec:=Prec);
		M_2 := SE_SmallPeriodMatrix(f,m:Prec:=Prec);
		
		
		// Compute the period matrices of y^m = f(x) to precision Prec digits (defined over complex field)

		C<I> := ComplexField(Prec);
		Cx<x> := PolynomialRing(C);
		f := x^6 + (6+28*I)*x^5 - 1;
		m := 10;
		Prec := 200;
		M_3 := SE_BigPeriodMatrix(f,m:Prec:=Prec);
		M_4 := SE_SmallPeriodMatrix(f,m:Prec:=Prec);
		

		// Compute the period matrices of y^m = f(x) to precision Prec digits defined via branch points

		Prec := 200;
		C<I> := ComplexField(Prec);
		Points := [ 5+3*I, 77, 1001, 2*I, -1-I ]; c := (1+I)^2;
		m := 4;
		M_5 := SE_BigPeriodMatrix(Points,m:Prec:=Prec,LeadingCoeff:=c);
		M_6 := SE_SmallPeriodMatrix(Points,m:Prec:=Prec,LeadingCoeff:=c);


		// Create a superelliptic curve object via univariate polynomial
		S := SE_Curve(f,m:Prec:=Prec);

		Optional parameters for intrinsic SE_Curve()
			Prec: 
				Prescribed precision for the complex superelliptic curve;
				Default: 20
			Small: 
				Compute small period matrix? (not needed for Abel-Jacobi);
				Default: true
			AbelJacobi: 
				Do precomputations for AbelJacobi map?;
				Default: true
			InftyPoints: 
				Compute S`AJM_InftyPoints (may take some time,see below); 
				Default: false
			IntegrationType: 
				Force use of double-exponential integration for m = 2; 
				Default: "Opt";
				Other: "DE","GC"
			

		// Create a superelliptic curve object via branch points and leading coefficient

		Points := [ 5+3*I, 77, 1001, 2*I, -1-I ]; c := (1+I)^2;
		S := SE_Curve(Points,m:Prec:=Prec,LeadingCoeff:=c);

		// Try some random superelliptic curve
		S := SE_RandomCurve(4,6:Prec:=100);

		// Extended printing
		Print(S:Extended);

		// Access computed data:
		S`DefiningPolynomial;
		S`Degree;
		S`Genus;
		S`Prec;
		S`ComplexField;
		S`RealField;
		S`HolomorphicDifferentials;
		S`BranchPoints;
		S`SpanningTree;
		S`IntegrationType;
		S`IntersectionMatrix;
		S`ElementaryIntegrals;
		S`BigPeriodMatrix;
		S`SmallPeriodMatrix;
		S`TreeMatrix;
		S`ReductionMatrix;
		S`Basepoint;			// Basepoint for Abel-Jacobi map: P_0 = <S`Branchpoints[1],0>
		S`AJM_RamificationPoints;	// Abel-Jacobi map of D = [P_k - P_0] for k = 1,..,deg(f)
		S`AJM_SumOfInftyPoints;		// Abel-Jacobi map of D = [ \sum P_{\infty} - gcd(m,deg(f)) P_0 ]
		S`AJM_InftyPoints;		// Abel-Jacobi map of D = [ P^{(k)}_{\infty} - P_0 ], k = 1,..,gcd(m,deg(f))

		// S`AJM_InftyPoints are not computed by default (depending on the curve this may take a while)
		// Compute Abel-Jacobi map of D = [ P^{(k)}_{\infty} - P_0 ]

		for k in [1..S`Degree[3]] do
			SE_AJM_InftyPoints(k,S);
		end for;


		// Abel-Jacobi map of divisors with image in \R^{2g} / Z^{2g}

		f := BernoulliPolynomial(10);
		m := 5;
		S := SE_Curve(f,m:Prec:=100); // Genus 16
			

		// Compute Abel-Jacobi map between P_1 and P_2

		P_1 := SE_RandomPoint(S);
		P_2 := SE_RandomPoint(S);
		A := SE_AbelJacobi(P_1,P_2,S);


		// Define a divisor via points on the curve
		
		C<I> := S`ComplexField;	
		x_1 := One(C); Fx_1 := SE_Branches(x_1,S);
		x_2 := 1+I;    Fx_2 := SE_Branches(x_2,S);
		D := SE_Divisor([<[x_1,Fx_1[1]],5>,<[x_1,Fx_1[3]],-3>,<[x_2,Fx_2[2]],10>],S);
		

		// Compute Abel-Jacobi map of D from Q_0

		Q_0 := [x_1,Fx_1[3]];
		A := SE_AbelJacobi(D,Q_0,S);
		

		// Diviors can also include infinite points

		D := SE_Divisor([<[1],1>,<[2],1>,<[3],1>,<[4],1>,<[5],1>],S);
		A := SE_AbelJacobi(D,Q_0,S);


		// Abel-Jacobi map of degree zero divisors indepently of a basepoint

		D := SE_Divisor([<P_1,-10>,<[4],10>],S);
		A := SE_AbelJacobi(D,S);


		// Or pick a random zero divisor on k points

		k := 5;
		D := SE_RandomDivisor(k,S:Zero:=true);
		A := SE_AbelJacobi(D,S);

