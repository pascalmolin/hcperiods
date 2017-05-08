/******************************************************************************

 Copyright (C) 2017 Christian Neurohr

 ******************************************************************************/


function SE_TestFamilies( m, n : Prec:=20 )
// Test superelliptic period matrices for different families of polynomials

	Small := true;
	AbelJacobi := true;
	InftyPoints := true;

	Qx<x> := PolynomialRing(Rationals());

	// Cyclotomic polynomials
	print "Testing cyclotomic polynomials...";
	C_n := x^n + 1;
	print "C_n:",C_n;
	time S := SE_Curve(C_n,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,InftyPoints:=InftyPoints);
	print "Check.";


	// Exponential series
	print "Testing exponential series...";
	E_n := &+[ x^k/Factorial(k) : k in [0..n] ];
	E_n_R := Qx!Reverse(Coefficients(E_n));
	print "E_n:",E_n;
	time S := SE_Curve(E_n,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,InftyPoints:=InftyPoints);
	print "E_n_R:",E_n_R;
	time S := SE_Curve(E_n_R,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,InftyPoints:=InftyPoints);
	print "Check.";	

	// Bernoulli polynomial
	print "Testing Bernoulli polynomials...";
	B_n := BernoulliPolynomial(n);
	B_n_R := Qx!Reverse(Coefficients(B_n));
	print "B_n:",B_n;
	time S := SE_Curve(B_n,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,InftyPoints:=InftyPoints);
	if Degree(B_n_R) ge 3 then
		print "B_n_R:",B_n_R;
		time S := SE_Curve(B_n_R,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,InftyPoints:=InftyPoints);
	end if;
	print "Check.";
	
	//  Legendre polynomials
	print "Testing Legendre polynomials...";
	P_n := (1/2^n) * &+[ Binomial(n,k)^2 * (x-1)^(n-k) * (x+1)^k : k in [0..n] ];
	P_n_R := Qx!Reverse(Coefficients(P_n));
	print "P_n:",P_n;
	time S := SE_Curve(P_n,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,InftyPoints:=InftyPoints);
	if Degree(P_n_R) ge 3 then
		print "P_n_R:",P_n_R;
		time S := SE_Curve(P_n_R,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,InftyPoints:=InftyPoints);
	end if; 
	print "Check.";
	
	// Laguerre polynomials
	print "Testing Laguerre polynomials...";
	L_n := &+[ Binomial(n,k) * (-1)^k / Factorial(k) * x^k : k in [0..n] ];
	L_n_R := Qx!Reverse(Coefficients(L_n));
	print "L_n:",L_n;
	time S := SE_Curve(L_n,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,InftyPoints:=InftyPoints);
	print "L_n_R:",L_n_R;
	time S := SE_Curve(L_n_R,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,InftyPoints:=InftyPoints);
	print "Check.";

	// Chebyshev polynomials
	print "Testing Chebyshev polynomials...";
	T_n := &+[ Binomial(n,2*k) * (x^2-1)^k * x^(n-2*k) : k in [0..Floor(n/2)] ];
	T_n_R := Qx!Reverse(Coefficients(T_n));
	print "T_n:",T_n;
	time S := SE_Curve(T_n,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,InftyPoints:=InftyPoints);
	if Degree(T_n_R) ge 3 then
		print "T_n_R:",T_n_R;
		time S := SE_Curve(T_n_R,m:Prec:=Prec,Small:=Small,AbelJacobi:=AbelJacobi,InftyPoints:=InftyPoints);
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


function SE_AJM_LongTest( m, n : Prec:=20, t:=10 )
// Long time test for Abel-Jacobi map
	//Err := SEC`Error;
	Err := 10^-Floor(Prec/2);
	for j in [1..t] do
		SEC := SE_RandomCurve(m,n:Prec:=Prec);
		for k in [1..t] do
			V := SE_AJM_Test_2(SEC);
			for l in [1..2*SEC`Genus] do
				if V[l][1] gt Err then
					print "SEC:",SEC;
					error Error("Abel-Jacobi map test 2 failed.");
				end if;
			end for;
		end for;
	end for;
	return true;
end function;


function SE_AJM_LongTest_Random( m, n : Prec:=20, t:=10 )
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















