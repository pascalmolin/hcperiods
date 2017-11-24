/******************************************************************************

 Copyright (C) 2017 Christian Neurohr

 ******************************************************************************/


function SE_TestFamilies( m, n : Prec := 30, Exact := true, Small := true, AbelJacobi := true, InftyPoints := false, IntegrationType := "Opt" )
// Test superelliptic period matrices for different families of polynomials

	if Exact then
		Kx<x> := PolynomialRing(Rationals());
	else
		C<I> := ComplexField(Prec);
		Kx<x> := PolynomialRing(C);
	end if;		

	// Cyclotomic polynomials
	print "Testing cyclotomic polynomials...";
	C_n := x^n + 1;
	print "C_n:",C_n;
	time S := SE_Curve(C_n,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,InftyPoints:=InftyPoints,IntegrationType:=IntegrationType);
	print "Check.";

	// Exponential series
	print "Testing exponential series...";
	E_n := &+[ x^k/Factorial(k) : k in [0..n] ];
	E_n_R := Kx!Reverse(Coefficients(E_n));
	print "E_n:",E_n;
	time S := SE_Curve(E_n,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,InftyPoints:=InftyPoints,IntegrationType:=IntegrationType);
	print "E_n_R:",E_n_R;
	time S := SE_Curve(E_n_R,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,InftyPoints:=InftyPoints,IntegrationType:=IntegrationType);
	print "Check.";	

	// Bernoulli polynomial
	print "Testing Bernoulli polynomials...";
	B_n := Kx!BernoulliPolynomial(n);
	B_n_R := Kx!Reverse(Coefficients(B_n));
	print "B_n:",B_n;
	time S := SE_Curve(B_n,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,InftyPoints:=InftyPoints,IntegrationType:=IntegrationType);
	if Degree(B_n_R) ge 3 then
		print "B_n_R:",B_n_R;
		time S := SE_Curve(B_n_R,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,InftyPoints:=InftyPoints,IntegrationType:=IntegrationType);
	end if;
	print "Check.";
	
	//  Legendre polynomials
	print "Testing Legendre polynomials...";
	P_n := (1/2^n) * &+[ Binomial(n,k)^2 * (x-1)^(n-k) * (x+1)^k : k in [0..n] ];
	P_n_R := Kx!Reverse(Coefficients(P_n));
	print "P_n:",P_n;
	time S := SE_Curve(P_n,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,InftyPoints:=InftyPoints,IntegrationType:=IntegrationType);
	if Degree(P_n_R) ge 3 then
		print "P_n_R:",P_n_R;
		time S := SE_Curve(P_n_R,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,InftyPoints:=InftyPoints,IntegrationType:=IntegrationType);
	end if; 
	print "Check.";
	
	// Laguerre polynomials
	print "Testing Laguerre polynomials...";
	L_n := &+[ Binomial(n,k) * (-1)^k / Factorial(k) * x^k : k in [0..n] ];
	L_n_R := Kx!Reverse(Coefficients(L_n));
	print "L_n:",L_n;
	time S := SE_Curve(L_n,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,InftyPoints:=InftyPoints,IntegrationType:=IntegrationType);
	print "L_n_R:",L_n_R;
	time S := SE_Curve(L_n_R,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,InftyPoints:=InftyPoints,IntegrationType:=IntegrationType);
	print "Check.";

	// Chebyshev polynomials
	print "Testing Chebyshev polynomials...";
	T_n := &+[ Binomial(n,2*k) * (x^2-1)^k * x^(n-2*k) : k in [0..Floor(n/2)] ];
	T_n_R := Kx!Reverse(Coefficients(T_n));
	print "T_n:",T_n;
	time S := SE_Curve(T_n,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,InftyPoints:=InftyPoints,IntegrationType:=IntegrationType);
	if Degree(T_n_R) ge 3 then
		print "T_n_R:",T_n_R;
		time S := SE_Curve(T_n_R,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,InftyPoints:=InftyPoints,IntegrationType:=IntegrationType);
	end if;
	print "Check.";

	return true;
end function;


function SE_AJM_Test_1( SEC : Ht := 10^2 )
// Weak test for Abel-Jacobi map of superelliptic curves
	m := SEC`Degree[1];
	n := SEC`Degree[2];
	C<I> := ComplexField((Ceiling((m/10))+1)*SEC`Prec+2*n);
	x1 := C!Random([-Ht..Ht])/Random([1..Ht]) + I*C!Random([-Ht..Ht])/Random([1..Ht]);
	x2 := C!Random([-Ht..Ht])/Random([1..Ht]) + I*C!Random([-Ht..Ht])/Random([1..Ht]);
	F1 := SE_Branches(x1,SEC);
	F2 := SE_Branches(x2,SEC);
	Points := [];
	for j in [1..m] do
		Append(~Points,<[x1,F1[j]],1>);
		Append(~Points,<[x2,F2[j]],-1>);;
	end for;
	D := SE_Divisor(Points,SEC:Check:=true);
	return SE_AbelJacobi(D,SEC);
end function;


function SE_AJM_Test_2( SEC : Ht := 10^2 )
// Strong test for Abel-Jacobi map of superelliptic curves 
	m := SEC`Degree[1];
	n := SEC`Degree[2];
	C<I> := ComplexField((Ceiling((m/10))+1)*SEC`Prec+2*n);
	y1 := C!Random([-Ht..Ht])/Random([1..Ht]) + I*C!Random([-Ht..Ht])/Random([1..Ht]);
	y2 := C!Random([-Ht..Ht])/Random([1..Ht]) + I*C!Random([-Ht..Ht])/Random([1..Ht]);
	g1 := SEC`ComplexPolynomial - y1^m;
	g2 := SEC`ComplexPolynomial - y2^m;
	F1 := Roots(g1);
	F2 := Roots(g2);
	Points := [];
	for j in [1..n] do
		Append(~Points,<[F1[j][1],y1],1>);
		Append(~Points,<[F2[j][1],y2],-1>);
	end for;
	D := SE_Divisor(Points,SEC:Check:=true);
	return SE_AbelJacobi(D,SEC);
end function;


function SE_AJM_LongTest( m, n : Prec:=40, t:=10 )
// Long time test for Abel-Jacobi map
	//Err := SEC`Error;
	Err := 10^-Floor(Prec/2);
	for j in [1..t] do
		SEC := SE_RandomCurve(m,n:Prec:=Prec);
		for k in [1..t] do
			V := SE_AJM_Test_2(SEC);
			for l in [1..2*SEC`Genus] do
				if V[l][1] gt Err then
					print "V[l][1]:",V[l][1];
					print "Err:",Real(Err);
					print "SEC:",SEC;
					error Error("Abel-Jacobi map test 2 failed.");
				end if;
			end for;
		end for;
	end for;
	return true;
end function;


function SE_AJM_LongTest_Random( m, n : Prec:=40, t:=10 )
// Test Abel-Jacobi map for random divisors
	for j in [1..t] do
		SEC := SE_RandomCurve(m,n:Prec:=Prec);
		for k in [1..t] do
			D := SE_RandomDivisor(t,SEC);
			P := SE_RandomPoint(SEC);
			V := SE_AbelJacobi(D,P,SEC);
		end for;
	end for;
	return true;
end function;

/* MAGMA error 
Kx<x> := PolynomialRing(Rationals());
B_n := Kx!BernoulliPolynomial(30);
B_n_R := Kx!Reverse(Coefficients(B_n));
A := AnalyticJacobian(ChangeRing(B_n_R,ComplexField(40)));
*/

intrinsic SE_TestIntegration(m::RngIntElt,n::RngIntElt:Prec:=30,Test:=false,Small:=false) -> BoolElt
{ -- }
	Kx<x> := PolynomialRing(Rationals());
	B_n := Kx!BernoulliPolynomial(n);
	B_n_R := Kx!Reverse(Coefficients(B_n));	

	
	if false then
	Timings_B_GJ := [];
	Timings_Br_GJ := [];
	Timings_B_DE := [];
	Timings_Br_DE := [];
	Timings_B_MAGMA := [];
	Timings_Br_MAGMA := [];
	P := [50,200,500,1000,2000];
	Q := [8,12,24,48,96];
	for p in P do
		Prec := p;
		C := ComplexField(Prec);
		for q in Q do
			B_n := Kx!BernoulliPolynomial(q);
			B_n_R := Kx!Reverse(Coefficients(B_n));
			print "p:",p;
			print "q:",q;
			"B_n,GJ:";
			t := Cputime();
			time S11 := SE_Curve(B_n,m:Prec:=Prec,Small:=Small,AbelJacobi:=false,InftyPoints:=false,IntegrationType:="GJ");
			Append(~Timings_B_GJ,[p,q,Cputime()-t]);
			"B_n,DE:";
			t := Cputime();
			time S12 := SE_Curve(B_n,m:Prec:=Prec,Small:=Small,AbelJacobi:=false,InftyPoints:=false,IntegrationType:="DE");
			Append(~Timings_B_DE,[p,q,Cputime()-t]);
			"Check isomorphism:";
			//bool1 := IsIsomorphicSmallPeriodMatrices(S11`SmallPeriodMatrix,S12`SmallPeriodMatrix); print "bool1:",bool1;
			"B_n_R,GJ:";
			t := Cputime();
			time S21 := SE_Curve(B_n_R,m:Prec:=Prec,Small:=Small,AbelJacobi:=false,InftyPoints:=false,IntegrationType:="GJ");
			Append(~Timings_Br_GJ,[p,q,Cputime()-t]);
			"B_n_R,DE:";
			t := Cputime();
			time S22 := SE_Curve(B_n_R,m:Prec:=Prec,Small:=Small,AbelJacobi:=false,InftyPoints:=false,IntegrationType:="DE");
			Append(~Timings_Br_DE,[p,q,Cputime()-t]);
			//bool2 := IsIsomorphicSmallPeriodMatrices(S21`SmallPeriodMatrix,S22`SmallPeriodMatrix); print "bool2:",bool2;
			if q notin [24,48,96] and p ne 2000 then
				"B_n,Magma:";
				t := Cputime();
				time A1 := AnalyticJacobian(ChangeRing(B_n,C));
				Append(~Timings_B_MAGMA,[p,q,Cputime()-t]);
				//bool3 := IsIsomorphicSmallPeriodMatrices(ChangeRing(S11`SmallPeriodMatrix,C),SmallPeriodMatrix(A1)); print "bool3:",bool3;
				"B_n_R,Magma:";
				t := Cputime();
				time A2 := AnalyticJacobian(ChangeRing(B_n_R,C));
				Append(~Timings_Br_MAGMA,[p,q,Cputime()-t]);
				//bool4 := IsIsomorphicSmallPeriodMatrices(ChangeRing(S21`SmallPeriodMatrix,C),SmallPeriodMatrix(A2)); print "bool4:",bool4;
			end if;
		end for;
	end for;
	//Timings := [Timings_B_GJ,Timings_B_DE,Timings_Br_GJ,Timings_Br_DE];
	Timings := [Timings_B_GJ,Timings_B_DE,Timings_Br_GJ,Timings_Br_DE,Timings_B_MAGMA,Timings_Br_MAGMA];
	return Timings;
	end if;

	if false then
	Timings_B_GJ := [];
	Timings_Br_GJ := [];
	Timings_B_DE := [];
	Timings_Br_DE := [];
	P := [50,200,500,1000,2000]; P := [2000];
	//Q := [[3,4,0],[4,6,1],[7,7,0],[10,8,1],[8,12,0],[24,8,0],[18,12,1],[17,24,0],[24,48,1],[18,96,0]];
	Q := [[24,48,1]];
	for p in P do
		Prec := p;
		C := ComplexField(Prec);
		for q in Q do
			B_n := Kx!BernoulliPolynomial(q[2]);
			B_n_R := Kx!Reverse(Coefficients(B_n));
			"###########################";
			print "p:",p;
			print "q:",q;
			m := q[1];
			if q[3] eq 0 then
				B := B_n;
			else
				B := B_n_R;
			end if;
			if true then
				"DE";
				t := Cputime();
				time S12 := SE_Curve(B,m:Prec:=Prec,Small:=Small,AbelJacobi:=false,InftyPoints:=false,IntegrationType:="DE");
				Append(~Timings_B_DE,<p,q,Cputime()-t>);
			end if;
			if false then
				bool1 := IsIsomorphicSmallPeriodMatrices(S11`SmallPeriodMatrix,S12`SmallPeriodMatrix); print "bool1:",bool1;
			end if;
			if true then
				"GJ:";
				t := Cputime();
				time S11 := SE_Curve(B,m:Prec:=Prec,Small:=Small,AbelJacobi:=false,InftyPoints:=false,IntegrationType:="GJ");
				Append(~Timings_B_GJ,<p,q,Cputime()-t>);
			end if;
		end for;
	end for;
	Timings := [Timings_B_GJ,Timings_B_DE];
	return Timings;
	end if;


	if false then
	Timings_B_GJ := [];
	Timings_Br_GJ := [];
	Timings_B_DE := [];
	Timings_Br_DE := [];
	Timings_B_MAGMA := [];
	Timings_Br_MAGMA := [];
	P := [50,200,500,1000,2000];
	Q := [8,12,24,48,96];
	for p in P do
		Prec := p;
		C := ComplexField(Prec);
		for q in Q do
			B_n := Kx!BernoulliPolynomial(q);
			B_n_R := Kx!Reverse(Coefficients(B_n));
			print "p:",p;
			print "q:",q;
			"B_n,GJ:";
			t := Cputime();
			time S11 := SE_Curve(B_n,m:Prec:=Prec,Small:=Small,AbelJacobi:=false,InftyPoints:=false,IntegrationType:="GJ");
			Append(~Timings_B_GJ,[p,q,Cputime()-t]);
			"B_n,DE:";
			t := Cputime();
			time S12 := SE_Curve(B_n,m:Prec:=Prec,Small:=Small,AbelJacobi:=false,InftyPoints:=false,IntegrationType:="DE");
			Append(~Timings_B_DE,[p,q,Cputime()-t]);
			//"Check isomorphism:";
			//bool1 := IsIsomorphicSmallPeriodMatrices(S11`SmallPeriodMatrix,S12`SmallPeriodMatrix); print "bool1:",bool1;
			"B_n_R,GJ:";
			t := Cputime();
			time S21 := SE_Curve(B_n_R,m:Prec:=Prec,Small:=Small,AbelJacobi:=false,InftyPoints:=false,IntegrationType:="GJ");
			Append(~Timings_Br_GJ,[p,q,Cputime()-t]);
			"B_n_R,DE:";
			t := Cputime();
			time S22 := SE_Curve(B_n_R,m:Prec:=Prec,Small:=Small,AbelJacobi:=false,InftyPoints:=false,IntegrationType:="DE");
			Append(~Timings_Br_DE,[p,q,Cputime()-t]);
		end for;
	end for;
	Timings := [Timings_B_GJ,Timings_B_DE,Timings_Br_GJ,Timings_Br_DE];
	//Timings := [Timings_B_GJ,Timings_B_DE,Timings_Br_GJ,Timings_Br_DE,Timings_B_MAGMA,Timings_Br_MAGMA];
	return Timings;
	end if;

	/*if Test then
		P := [200];
		Lxy<x,y> := PolynomialRing(Rationals(),2);
		g := y^m - Evaluate(B_n,x);
		gr := y^m - Evaluate(B_n_R,x);
		PM1 := RS_PeriodMatrix(g:Prec:=Prec);
		PM2 := RS_PeriodMatrix(gr:Prec:=Prec);
		time A1 := AnalyticJacobian(ChangeRing(B_n,C));
		time A2 := AnalyticJacobian(ChangeRing(B_n_R,C));
		bool1 := IsIsomorphicSmallPeriodMatrices(SmallPeriodMatrix(A1),SE_SmallPeriodMatrix(B_n,m:Prec:=Prec));
		bool2 := IsIsomorphicSmallPeriodMatrices(SmallPeriodMatrix(A2),SE_SmallPeriodMatrix(B_n_R,m:Prec:=Prec));
		bool3 := IsIsomorphicSmallPeriodMatrices(SmallPeriodMatrix(A1),PM1);
		bool4 := IsIsomorphicSmallPeriodMatrices(SmallPeriodMatrix(A2),PM2);
		print "b1",bool1;
		print "b2",bool2;
		print "b3",bool3;
		print "b4",bool4;
		//bool3 := IsIsomorphicBigPeriodMatrices(BigPeriodMatrix(A1),S11`BigPeriodMatrix);
		bool5 := IsIsomorphicBigPeriodMatrices(S11`BigPeriodMatrix,S21`BigPeriodMatrix);
		//bool5 := IsIsomorphicBigPeriodMatrices(BigPeriodMatrix(A2),S12`BigPeriodMatrix);
		bool6 := IsIsomorphicBigPeriodMatrices(S12`BigPeriodMatrix,S22`BigPeriodMatrix);
		//bool7 := IsIsomorphicBigPeriodMatrices(S11`BigPeriodMatrix,RS_PeriodMatrix(g:Small:=false,Prec:=Prec+2)); print "b7",bool7;
		//bool8 := IsIsomorphicBigPeriodMatrices(S12`BigPeriodMatrix,RS_PeriodMatrix(gr:Small:=false,Prec:=Prec+2)); print "b8",bool8;
		print "b5",bool5;
		print "b6",bool6;			
	end if;*/

	if true then
	Timings_GJ := [];
	Timings_DE := [];
	P := [50,200,500,1000,2000];
	Q := [[2,8],[5,5],[4,6],[7,7],[6,8],[8,12]];
	for p in P do
		Prec := p;
		C := ComplexField(Prec);
		for q in Q do
			m := q[1];
			n := q[2];
			B_n := Kx!BernoulliPolynomial(n);
			print "p:",p;
			print "q:",q;
			S1 := SE_Curve(B_n,m:Prec:=Prec,Small:=false,AbelJacobi:=true,InftyPoints:=false,IntegrationType:="DE");
			S2 := SE_Curve(B_n,m:Prec:=Prec,Small:=false,AbelJacobi:=true,InftyPoints:=false,IntegrationType:="GJ");
			print "Computing Abel-Jacobi map to infinity:";
			"with double-exponential:";
			t := Cputime();
			for k in [1..S1`Degree[3]] do
				time SE_AJM_InftyPoints(k,S1);
			end for;
			Append(~Timings_DE,<p,q,Cputime()-t>);
			print "Total time(DE):",Cputime()-t;
			T1 := (m/Gcd(m,n))*&+S1`AJM_InftyPoints;
			assert &and[ v - Round(v) lt 10^-40 : v in Eltseq(T1) ];
			//print "TEST:",T1;		
			"with Gauss-Jacobi:";
			t := Cputime();
			for k in [1..S2`Degree[3]] do
				time SE_AJM_InftyPoints(k,S2);
			end for;
			Append(~Timings_GJ,<p,q,Cputime()-t>);
			print "Total time(GJ):",Cputime()-t;	
			T2 := (m/Gcd(m,n))*&+S2`AJM_InftyPoints;
			assert &and[ v - Round(v) lt 10^-40 : v in Eltseq(T2) ];
			//print "TEST:",T2;	
		end for;
	end for;
	Timings := [Timings_DE,Timings_GJ];
	return Timings;
	end if;
	
	
	/*t := 10;
	D := SE_RandomDivisor(t,S1);
	P := SE_RandomPoint(S1);
	print "Computing Abel-Jacobi map of divisor...";
	//print "and basepoint:",P;
	print "using double-exponential integration:";
	time V1 := SE_AbelJacobi(D,P,S1);
	print "using Gauss-Jacobi integration:";
	time V2 := SE_AbelJacobi(D,P,S2);*/
	return true;
end intrinsic;










