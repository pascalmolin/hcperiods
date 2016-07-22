//Import global settings
import "globalsettings.m": RS_Config;
import "superelliptic.m": RS_ImSgn, RS_NRootAC;
import "comparefunctions.m": RS_CompareFldComElt;


intrinsic RS_RootTest( Z::[FldComElt] ) -> RngIntElt
{ - }
	a := Z[1]; b := Z[2]; c := Z[3]; d := Z[4]; l := Random([3..10]); N := 3;
	assert Re(b-a) ge 0;
	assert Re(d-c) ge 0;
	Zeta_k := (-((1/(d-c))^l))^(1/N) / (-((1/(b-a))^l))^(1/N) / ((b-a)/(d-c))^(l/N);
	return Zeta_k;
end intrinsic;
	 


intrinsic RS_IntersectionTest( n::RngIntElt : N := 2 ) -> RngIntElt
{ - }
	C<i> := ComplexField(20);
	/*
	for l in [1..n] do
		k := 10^l;
		//k := n;
		Pts := [i,1+(1/k)+i,1-i,1+i+(1/k)*i,8+8*i,-7/i];
		E1 := <2,1>; E2 := <3,4>;
		s := RS_InnerIntersection(E1,E2,Pts,N);
		"############################################### neeeext";
	end for;
	*/

	a := -Random([1..n]) - i*Random([1..n]);
	b := -Random([1..n]) + i*Random([1..n]);
	c := Random([1..n]) + i*Random([1..n]);
	d := Random([1..n]) - i*Random([1..n]);

	Pts := [a,b,c,d];
	E1 := <1,3>; E2 := <2,4>;
	s := RS_InnerIntersection(E1,E2,Pts,N);
	print "s:",s;
	E3 := <3,1>; E4 := <2,4>;
	s := RS_InnerIntersection(E3,E4,Pts,N);
	print "s:",s;
	return 0;
end intrinsic;



intrinsic RS_CircleTest( X::RieSrf : r := 1/10, j := 1 ) -> RngIntElt
{ - }
	RS_PeriodMatrix(X); f := X`DefiningPolynomial; N := Degree(f,2); k := 1; g := X`Genus;
	Pts := Sort(X`DiscriminantPoints,RS_CompareFldComElt); print "Pts:",Pts; print "j:",j;
	a := Pts[j]; print "a:",a;
	b := Pts[j+1]; print "b:",b;
	ls0 := RS_LineSegment(a,b);
	fc1 := RS_FullCircle(a-r,a:Orientation:=1);
	fc2 := RS_FullCircle(b-r,b:Orientation:=1);
	Okay1, IntPts1 := RS_IntersectionPoints( ls0, fc1 ); print "IntPts1:",IntPts1;
 	Okay2, IntPts2 := RS_IntersectionPoints( ls0, fc2 ); print "IntPts2:",IntPts2;
	assert &and[Okay1,Okay2]; //assert &and[#IntPts1 eq 1,#IntPts2 eq 1];
	FC1 := RS_Chain(RS_FullCircle(IntPts1[1],a:Orientation:=1)); RS_ChainIntegrals(X,FC1); print "FC1:",FC1;
	FC2 := RS_Chain(RS_FullCircle(IntPts2[1],b:Orientation:=1)); RS_ChainIntegrals(X,FC2);print "FC2:",FC2;
	LS := RS_Chain(RS_LineSegment(FC1`StartPt,FC2`StartPt)); RS_ChainIntegrals(X,LS); print "LS:",LS;
	LSInv := LS^(-1);
	Chain := LS * FC2^(N-1) * LSInv * FC1;
	//Chain := LS * FC2 * LSInv * FC1;
	//Chain := LS * FC2^k * LS^(-1) * FC1^(N-k);
	RS_ChainIntegrals(X,Chain);
	print "FC1`Integrals:",FC1`Integrals;
	print "FC2`Integrals:",FC2`Integrals;
	print "LS`Integrals:",LS`Integrals;
	print "LSInv`Integrals:",LSInv`Integrals;
	SE_DFF := RS_SuperellipticDifferentials(f);
	SE_Integrals := RS_SEIntegration( f, ls0, Pts, SE_DFF );
	print "SE_Integrals:",SE_Integrals;
	print "Chain`Integrals:",Chain`Integrals;
	print "Chain`Permutation:",Chain`Permutation;
	Fails := [];
	for j in [1..g] do
		thisw := false;
		for k in [1..g] do
			Okay, Sigma := RS_Permutation(SE_Integrals[j],Chain`Integrals[k]`Entries,10^-20:Distinct:=false);
			if Okay then
				print "Sigma:",Sigma;
				print "k:",k;
				thisw := true;
				continue;
			end if;		
		end for;
		
		if thisw eq false then
			print "j:",j;
			Append(~Fails,SE_DFF[j]);
		end if;
	end for;
	return Fails;
end intrinsic;	


intrinsic RS_TestLimits( Tau::FldReElt : N := 2 ) -> RngIntElt
{ - }
	Pi := RS_GetGlobalPi();
	C<i> := RS_GetGlobalComplexField_20();
	assert &and[-Pi/2 lt Tau, Tau lt Pi/2];
	print "Tau:",Tau;
	x := 0;
	for j in [1..10] do
		x +:= j;
		print "x:",x;
		print "Arg(1-Tanh(x+i*Tau)):",Arg((1-Tanh(x+i*Tau))^(2/N));
		print "Arg(1-Tanh(x+i*Tau)):",Arg((1-Tanh(x-i*Tau))^(2/N));
		//print "Arg(1+Tanh(x+i*Tau)):",Arg((1+Tanh(x+i*Tau))^(2/N));
		//print "Arg(1+Tanh(x-i*Tau)):",Arg((1+Tanh(x-i*Tau))^(2/N));
		print "Arg(Cosh(x+i*Tau)):", Arg((Cosh(x+i*Tau))^(2/N));
		print "Arg(Cosh(x-i*Tau)):", Arg((Cosh(x-i*Tau))^(2/N));
		print "Arg(Cosh(-x+i*Tau)):", Arg((Cosh(-x+i*Tau))^(2/N));
		print "Arg(Cosh(-x-i*Tau)):", Arg((Cosh(-x-i*Tau))^(2/N));
		"###########################";
	end for;
	return 0;
end intrinsic;


intrinsic RS_TestRoot1( z::FldComElt : N:= 2 ) -> RngIntElt
{ - }
	assert Im(z) ne 0;
	print "z:",z;
	print "N:",N;
	z1 := z^(1/N);
	r,phi := ComplexToPolar(z);
	print "r:",r;
	print "phi:",phi;
	"######'";
	print "z^(1/N):",z1;
	z2 := r^(1/N) * ( Cos(phi/N) + Parent(z).1 * Sin(phi/N) );
	print "r^(1/N) * ( Cos(phi/N) + i * Sin(phi/N) ):",r^(1/N) * ( Cos(phi/N) + Parent(z).1 * Sin(phi/N) );
	return Abs(z1-z2) le RS_GetGlobalError();
end intrinsic;



intrinsic RS_TestSqrt( z1::FldComElt, z2::FldComElt : N := 2 ) -> RngIntElt
{ - }
	print "####";
	print "z1:",z1; 
	print "z2:",z2;
	print "Arg(z1):",Arg(z1);
	print "Arg(z2):",Arg(z2);
	print "z1*z2:",z1*z2; 
	print "(z1*z2)^(1/N):",(z1*z2)^(1/N);
	print "Arg(z1*z2):",Arg(z1*z2);
	print "Arg((z1*z2)^(1/N)):",Arg((z1*z2)^(1/N));
	print "z1^(1/N) * z2^(1/N):",z1^(1/N) * z2^(1/N);
	print "Arg(z1^(1/N) * z2^(1/N)):",Arg(z1^(1/N) * z2^(1/N));  
	print "####";
	return z1^(1/N) * z2^(1/N)/(z1*z2)^(1/N);
end intrinsic;

intrinsic RS_TestRoot2( z::FldComElt, N::RngIntElt ) -> RngIntElt
{ - }
	assert N ge 2; C_1 := One(Parent(z)); Pi := RS_GetGlobalPi(); i := Parent(z).1;
	if N eq 2 then
		if Im(z) eq 0 and Re(z) lt 0 then
			Fact := -i;
		else
			Fact := 1;
		end if;
	else
		if RS_ImSgn(z) gt 0 then
			Fact := 1;
		else
			Fact := 1/Exp(2*Pi*Parent(z).1);
		end if;
	end if;
	print "z:",z;
	print "Fact:",Fact;
	print "(-z)^(1/N):",(-z)^(1/N);
	print "Fact * (-z)^(1/N):",Fact * (-z)^(1/N);
	print "z^(1/N):",z^(1/N);
	return 0;
end intrinsic;

intrinsic RS_TestConstants( t::FldReElt, Tau::FldReElt ) -> RngIntElt
{ - }
	Pi := RS_GetGlobalPi();
	C<i> := RS_GetGlobalComplexField_20();
	//assert &and[0 lt Tau, Tau lt Pi/2];
	//LHS := (1 + Tanh(t - i*Tau)) * ( 1 - Tanh(t - i*Tau)); print "LHS:",LHS;
	//RHS := (1 / (Cosh(t - i*Tau))^2); print "RHS:",RHS;
	N := Random([2..100]); print "N:",N; 
	j := Random([1..N-1]); print "j:",j;
	c1 := Exp(-Pi*i/N); print "c1:",c1;
	c2 := Exp(Pi*i/N); print "c2:",c2;
	LHS1 := 1 / ((Tanh(t + i*Tau)+1)^(j/N) * (Tanh(t+i*Tau)-1)^(j/N)); print "LHS1:",LHS1;
	RHS1 := c1^j * Cosh(t + i*Tau)^(2*j/N); print "RHS1:",RHS1;
	LHS2 := 1 / ((Tanh(t - i*Tau)+1)^(j/N) * (Tanh(t-i*Tau)-1)^(j/N)); print "LHS2:",LHS2;
	RHS2 := c2^j * Cosh(t - i*Tau)^(2*j/N); print "RHS2:",RHS2;
	
	"#####";
	//c1 := 1/((Tanh(t + i*Tau)+1)^(1/N) * (Tanh(t+i*Tau)-1)^(1/N) * Cosh(t+i*Tau)^(2/N)); print "c1:",c1; print "Abs(c1):",Abs(c1);
	//c2 := 1/((Tanh(t - i*Tau)+1)^(1/N) * (Tanh(t-i*Tau)-1)^(1/N) * Cosh(t-i*Tau)^(2/N)); print "c2:",c2; print "Abs(c2):",Abs(c2);
	return 0;
end intrinsic;


intrinsic RS_TestSpeedNRoot( d::RngIntElt, N::RngIntElt : Prec := 20, Height := 100 ) -> RngIntElt
{ Doesn't seem to work?!? }
	Z := RS_RandomPts(d:Prec:=Prec, Height:=Height); Zeta := RS_GetGlobalZeta();
	x := Random([-Height..Height])/Height;
	print "x:",x;
	print "Time (Long method):";
	time A1 := RS_NRootAC(x,Z,Zeta,N:long:=true);
	print "Time (Short method):";
	time A2 := RS_NRootAC(x,Z,Zeta,N:long:=false);
	print "A1:",A1;
	print "A2:",A2;
	print "Difference:",Abs(A1-A2);
	return 0;
end intrinsic;


intrinsic RS_RandomPts( n::RngIntElt : Prec:= 20, Height := 100 ) -> SeqEnum[FldComElt]
{ - }
	Z := [];
	C<i> := ComplexField(Prec);
	for j in [1..n] do
		z := Random([-Height..Height])/Height + Random([-Height..Height])/Height*i;
		Append(~Z,z);
	end for;
	return Z;
end intrinsic;



intrinsic RS_SEPMTest( s::RngMPolElt : Prec := 20 ) -> BoolElt, Matrix
{ Tests if period matrices of superelliptic curves are equivalent }
	print "Computing PM1 with general algorithm:";
	time PM1 := RS_PeriodMatrix(s:Prec:=Prec);
	print "Computing PM1 with special case algorithm:";
	time PM2 := RS_SEPM(s:Prec:=Prec);
	print "Searching symplectic transformation...";
	bool, S := IsIsomorphicSmallPeriodMatrices(PM1,PM2);
	return bool, S;
end intrinsic;


intrinsic RS_SEPMLTTest( N::RngIntElt, d::RngIntElt, T::RngIntElt: Ht := 100, Prec := 20 ) -> BoolElt
{ Tests if period matrices of superelliptic curves are equivalent }
	assert T gt 0;
	for t in [1..T] do
		s := RS_RandomSuperellipticCurve(N,d:Ht:=Ht);
		bool, S := RS_SEPMTest(s : Prec:= Prec);
		print "Test successful:",bool;
		print "Transformation:",S;
		if bool eq false then
			return s;
		end if;
	end for;
	return true;
end intrinsic;
























	



