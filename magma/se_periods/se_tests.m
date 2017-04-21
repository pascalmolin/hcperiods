/******************************************************************************

 Copyright (C) 2017 Christian Neurohr

 ******************************************************************************/


intrinsic SE_TestFamilies(m_max::RngIntElt,n_max::RngIntElt : Prec:=20 ) -> RngIntElt
{ Test superelliptic period matrices for different families of polynomial }
	
	assert &and[m_max ge 2, n_max ge 3];	

	Qx<x> := PolynomialRing(Rationals());

	// Cyclotomic polynomials
	print "Testing cyclotomic polynomials...";
	for n in [3..n_max] do
		C_n := x^n + 1;
		print "C_n:",C_n;
		for m in [2..m_max] do
			time S := SE_Curve(C_n,m:Prec:=Prec,Small:=true,InfinitePoints:=true);
		end for;
	end for;
	print "Check.";


	// Exponential series
	print "Testing exponential series...";
	for n in [3..n_max] do
		E_n := &+[ x^k/Factorial(k) : k in [0..n] ];
		print "E_n:",E_n;
		for m in [2..m_max] do
			time S := SE_Curve(E_n,m:Prec:=Prec,Small:=true,InfinitePoints:=true);
		end for;
	end for;
	print "Check.";	

	// Bernoulli polynomial
	print "Testing Bernoulli polynomials...";
	for n in [3..n_max] do
		B_n := BernoulliPolynomial(n);
		print "B_n:",B_n;
		for m in [2..m_max] do
			time S := SE_Curve(B_n,m:Prec:=Prec,Small:=true,InfinitePoints:=true);
		end for;
	end for;
	print "Check.";
	
	//  Legendre polynomials
	print "Testing Legendre polynomials...";
	for n in [3..n_max] do
		P_n := (1/2^n) * &+[ Binomial(n,k)^2 * (x-1)^(n-k) * (x+1)^k : k in [0..n] ];
		print "P_n:",P_n;
		for m in [2..m_max] do
			time S := SE_Curve(P_n,m:Prec:=Prec,Small:=true,InfinitePoints:=true); 
		end for;
	end for;
	print "Check.";
	
	// Laguerre polynomials
	print "Testing Laguerre polynomials...";
	for n in [3..n_max] do
		L_n := &+[ Binomial(n,k) * (-1)^k / Factorial(k) * x^k : k in [0..n] ];
		print "L_n:",L_n;
		for m in [2..m_max] do	
			time S := SE_Curve(L_n,m:Prec:=Prec,Small:=true,InfinitePoints:=true);
		end for;
	end for;
	print "Check.";

	// Chebyshev polynomials
	print "Testing Chebyshev polynomials...";
	for n in [3..n_max] do
		T_n := &+[ Binomial(n,2*k) * (x^2-1)^k * x^(n-2*k) : k in [0..Floor(n/2)] ];
		print "T_n:",T_n;
		for m in [2..m_max] do
			time S := SE_Curve(T_n,m:Prec:=Prec,Small:=true,InfinitePoints:=true);
		end for;
	end for;
	print "Check.";

	return true;
end intrinsic;


intrinsic SE_AJM_Test_1( SEC::SECurve : Ht := 1000 ) -> RngIntElt
{ Weak test for Abel-Jacobi map of superelliptic curves }
	C<i> := ComplexField(2*Precision(SEC`ComplexField));
	x1 := C!Random([-Ht..Ht])/Random([1..Ht]) + i*C!Random([-Ht..Ht])/Random([1..Ht]);
	x2 := C!Random([-Ht..Ht])/Random([1..Ht]) + i*C!Random([-Ht..Ht])/Random([1..Ht]);
	F1 := SE_PrincipalBranch(SEC,x1:Fiber);
	F2 := SE_PrincipalBranch(SEC,x2:Fiber);
	Points := [];
	for j in [1..SEC`Degree[1]] do
		Append(~Points,<[x1,F1[j]],1>);
		Append(~Points,<[x2,F2[j]],-1>);;
	end for;
	D := SE_Divisor(Points,SEC:Check:=true);	
	return SE_AbelJacobi(D,SEC);
end intrinsic;





intrinsic SE_AJM_Test_2( SEC::SECurve : Ht := 1000 ) -> RngIntElt
{ Strong test for Abel-Jacobi map of superelliptic curves }
	C<i> := ComplexField(2*Precision(SEC`ComplexField));
	m := SEC`Degree[1];
	y1 := C!Random([-Ht..Ht])/Random([1..Ht]) + i*C!Random([-Ht..Ht])/Random([1..Ht]); //print "y1:",y1;
	y2 := C!Random([-Ht..Ht])/Random([1..Ht]) + i*C!Random([-Ht..Ht])/Random([1..Ht]); //print "y2:",y2;
	g1 := SEC`ComplexPolynomial - y1^m;
	g2 := SEC`ComplexPolynomial - y2^m;
	F1 := Roots(g1);
	F2 := Roots(g2);
	Points := [];
	for j in [1..SEC`Degree[2]] do
		Append(~Points,<[F1[j][1],y1],1>);
		Append(~Points,<[F2[j][1],y2],-1>);
	end for;
	D := SE_Divisor(Points,SEC:Check:=true);
	return SE_AbelJacobi(D,SEC);
end intrinsic;


intrinsic SE_AJM_LongTest( m::RngIntElt, n::RngIntElt : Prec:=20, t := 10 ) -> RngIntElt
{ Long time test }
	for j in [1..t] do
		s := RS_RandomSuperellipticCurve(m,n);
		SEC := SE_Curve(s:Prec:=Prec);
		for k in [1..t] do
			V := SE_AJM_Test_2(SEC);
			print "V:",V;
			for l in [1..2*SEC`Genus] do
				//if V[l][1] gt 10^(-Prec+10) then
				if V[l][1] gt SEC`Error then
					print "SEC:",SEC;
					error Error("Abel-Jacobi map test 2 failed.");
				end if;
			end for;
		end for;
	end for;
	return true;
end intrinsic;
