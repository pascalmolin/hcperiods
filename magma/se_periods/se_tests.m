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
			time S := SE_Curve(C_n,m:Prec:=Prec);
		end for;
	end for;
	print "Check.";


	// Exponential series
	print "Testing exponential series...";
	for n in [3..n_max] do
		E_n := &+[ x^k/Factorial(k) : k in [0..n] ];
		print "E_n:",E_n;
		for m in [2..m_max] do
			time S := SE_Curve(E_n,m:Prec:=Prec);
		end for;
	end for;
	print "Check.";	

	// Bernoulli polynomial
	print "Testing Bernoulli polynomials...";
	for n in [3..n_max] do
		B_n := BernoulliPolynomial(n);
		print "B_n:",B_n;
		for m in [2..m_max] do
			time S := SE_Curve(B_n,m:Prec:=Prec);
		end for;
	end for;
	print "Check.";
	
	//  Legendre polynomials
	print "Testing Legendre polynomials...";
	for n in [3..n_max] do
		P_n := (1/2^n) * &+[ Binomial(n,k)^2 * (x-1)^(n-k) * (x+1)^k : k in [0..n] ];
		print "P_n:",P_n;
		for m in [2..m_max] do
			time S := SE_Curve(P_n,m:Prec:=Prec); 
		end for;
	end for;
	print "Check.";
	
	// Laguerre polynomials
	print "Testing Laguerre polynomials...";
	for n in [3..n_max] do
		L_n := &+[ Binomial(n,k) * (-1)^k / Factorial(k) * x^k : k in [0..n] ];
		print "L_n:",L_n;
		for m in [2..m_max] do	
			time S := SE_Curve(L_n,m:Prec:=Prec);
		end for;
	end for;
	print "Check.";

	// Chebyshev polynomials
	print "Testing Chebyshev polynomials...";
	for n in [3..n_max] do
		T_n := &+[ Binomial(n,2*k) * (x^2-1)^k * x^(n-2*k) : k in [0..Floor(n/2)] ];
		print "T_n:",T_n;
		for m in [2..m_max] do
			time S := SE_Curve(T_n,m:Prec:=Prec);
		end for;
	end for;
	print "Check.";

	return true;
end intrinsic;
