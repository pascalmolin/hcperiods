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

	Pts := [a,b,c,d,0,i,1];
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
	print "DFF:",X`BasisOfDifferentialsFirstKind;
	RS_PeriodMatrix(X); f := X`DefiningPolynomial; N := Degree(f,2); k := 1; g := X`Genus;
	Pts := Sort(X`DiscriminantPoints,RS_CompareFldComElt); print "Pts:",Pts; print "j:",j;
	a := Pts[j]; print "a:",a;
	b := Pts[j+1]; print "b:",b;
	assert Re(a) le Re(b);
	ls0 := RS_LineSegment(a,b);
	fc1 := RS_FullCircle(a-r,a:Orientation:=1);
	fc2 := RS_FullCircle(b-r,b:Orientation:=1);
	Okay1, IntPts1 := RS_IntersectionPoints( ls0, fc1 ); print "IntPts1:",IntPts1;
 	Okay2, IntPts2 := RS_IntersectionPoints( ls0, fc2 ); print "IntPts2:",IntPts2;
	assert &and[Okay1,Okay2]; //assert &and[#IntPts1 eq 1,#IntPts2 eq 1];

	// Define full circles
	FC1 := RS_Chain(RS_FullCircle(IntPts1[1],a:Orientation:=1)); 
	FC2 := RS_Chain(RS_FullCircle(IntPts2[1],b:Orientation:=1));
	RS_ChainIntegrals(X,FC1); 
	RS_ChainIntegrals(X,FC2);
	print "FC1:",FC1;
	print "FC2:",FC2;

	// Define line segments
	LS := RS_Chain(RS_LineSegment(FC1`StartPt,FC2`StartPt));
	RS_ChainIntegrals(X,LS);
	print "LS:",LS;
	LSInv := LS^(-1);

	// Define chain
	Chain := LS * FC2 * LSInv * FC1^-1;
	//Chain := LS * FC2 * LSInv * FC1^(N-1);
	//Chain := LS * FC2^-1 * LSInv * FC1;
	//Chain := LSInv * FC1 * LS * FC2^-1;
	


	//Chain := LSInv * FC1 * LS * FC2^(N-1);
	//Chain := LS * FC2^k * LS^(-1) * FC1^(N-k);
	//Chain := LS * FC2^(N-1) * LSInv * FC1;

	RS_ChainIntegrals(X,Chain);
	print "FC1`Integrals:",FC1`Integrals;
	print "FC2`Integrals:",FC2`Integrals;
	print "LS`Integrals:",LS`Integrals;
	print "LSInv`Integrals:",LSInv`Integrals;
	SE_DFF := RS_SEDifferentials(f);
	DFF := &cat[ DFF_i : DFF_i in SE_DFF ];
	SE_Integrals := RS_SEIntegration( f, ls0, Pts, DFF );
	print "LS`Perm:",LS`Permutation;
	print "FC1Perm:",FC1`Permutation;
	print "FC2Perm:",FC2`Permutation;
	print "SE_Integrals:",SE_Integrals;
	print "Chain`Integrals:",Chain`Integrals;
	print "Chain`Permutation:",Chain`Permutation;
	Fails := [];
	for j in [1..g] do
		thisw := false;
		for k in [1..g] do
			Okay, Sigma := RS_Permutation(SE_Integrals[j],Chain`Integrals[k]`Entries,10^-5);
			if Okay then
				print "Sigma:",Sigma;
				print "k:",k;
				thisw := true;
				continue;
			end if;		
		end for;
		
		if thisw eq false then
			print "j:",j;
			Append(~Fails,DFF[j]);
		end if;
	end for;
	return Fails;
end intrinsic;	


intrinsic RS_TestLimits( Tau::FldReElt : N := 2 ) -> RngIntElt
{ - }
	Pi := RS_GetGlobalPi();
	C<i> := RS_GetGlobalComplexField_20();
	//assert &and[-Pi/2 lt Tau, Tau lt Pi/2];
	print "Tau:",Tau;
	x := 0;
	for j in [1..10] do
		x +:= 1;
		print "x:",x;
		print "Arg(1-Tanh(x+i*Tau)):",Arg((1-Tanh(x+i*Tau))^(2/N));
		print "Arg(1-Tanh(x-i*Tau)):",Arg((1-Tanh(x-i*Tau))^(2/N));
		print "Arg(1-Tanh(-x+i*Tau)):",Arg((1-Tanh(-x+i*Tau))^(2/N));
		print "Arg(1-Tanh(-x-i*Tau)):",Arg((1-Tanh(-x-i*Tau))^(2/N));
		print "Arg(1+Tanh(x+i*Tau)):",Arg((1+Tanh(x+i*Tau))^(2/N));
		print "Arg(1+Tanh(x-i*Tau)):",Arg((1+Tanh(x-i*Tau))^(2/N));
		print "Arg(1+Tanh(-x-i*Tau)):",Arg((1+Tanh(-x-i*Tau))^(2/N));
		print "Arg(1+Tanh(-x+i*Tau)):",Arg((1+Tanh(-x+i*Tau))^(2/N));
		print "Arg(Cosh(x+i*Tau)):", Arg((Cosh(x+i*Tau))^(2/N));
		print "Arg(Cosh(x-i*Tau)):", Arg((Cosh(x-i*Tau))^(2/N));
		print "Arg(Cosh(-x+i*Tau)):", Arg((Cosh(-x+i*Tau))^(2/N));
		print "Arg(Cosh(-x-i*Tau)):", Arg((Cosh(-x-i*Tau))^(2/N));
		"###########################";
		print "(4/N)*Tau:",(4/N)*Tau;
		print "(2/N)*Tau:",(2/N)*Tau;
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
	LHS := (1 + Tanh(t - i*Tau)) * ( 1 - Tanh(t - i*Tau)); print "LHS:",LHS;
	RHS := (1 / (Cosh(t - i*Tau))^2); print "RHS:",RHS;
	LHS := (1 + Tanh(t + i*Tau)) * ( 1 - Tanh(t + i*Tau)); print "LHS:",LHS;
	RHS := (1 / (Cosh(t + i*Tau))^2); print "RHS:",RHS;
	N := Random([2..100]); print "N:",N; 
	j := Random([1..N-1]); print "j:",j;
	c1 := Exp(-Pi*i/N); print "c1:",c1;
	c2 := Exp(Pi*i/N); print "c2:",c2;
	LHS1 := 1 / ((Tanh(t + i*Tau)+1)^(j/N) * (Tanh(t+i*Tau)-1)^(j/N)); print "LHS1:",LHS1;
	RHS1 := c1^j * Cosh(t + i*Tau)^(2*j/N); print "RHS1:",RHS1;
	LHS2 := 1 / ((Tanh(t - i*Tau)+1)^(j/N) * (Tanh(t-i*Tau)-1)^(j/N)); print "LHS2:",LHS2;
	RHS2 := c2^j * Cosh(t - i*Tau)^(2*j/N); print "RHS2:",RHS2;
	/*
	print "Arg(Tanh(t+i*Tau)+1):",Arg(Tanh(t+i*Tau)+1);
	print "Arg(Tanh(t+i*Tau)-1):",Arg(Tanh(t+i*Tau)-1);
	print "Arg(Tanh(t-i*Tau)+1):",Arg(Tanh(t-i*Tau)+1);
	print "Arg(Tanh(t-i*Tau)-1):",Arg(Tanh(t-i*Tau)-1);
	print "Arg((Tanh(t-i*Tau)-1)*(Tanh(t-i*Tau)+1)):",Arg((Tanh(t-i*Tau)-1)*(Tanh(t-i*Tau)+1));
	print "Arg(1/((Tanh(t-i*Tau)-1)*(Tanh(t-i*Tau)+1))):",Arg(1/((Tanh(t-i*Tau)-1)*(Tanh(t-i*Tau)+1)));
	print "Arg(c1):",Arg(c1);
	print "Arg(c2):",Arg(c2);
	print "Arg(Cosh(t-i*Tau)):",Arg(Cosh(t-i*Tau));
	print "Arg(Cosh(t+i*Tau)):",Arg(Cosh(t+i*Tau));
	print "(Cosh(t - i*Tau))^2:",(Cosh(t - i*Tau))^2;
	print "-(Cosh(t - i*Tau))^2:",-(Cosh(t - i*Tau))^2;
	print "Arg((Cosh(t + i*Tau))^2):",Arg(Cosh(t + i*Tau)^2);
	print "Arg(-(Cosh(t + i*Tau))^2):",Arg(-Cosh(t + i*Tau)^2);
	*/
	"##########";
	//c1 := 1/((Tanh(t + i*Tau)+1)^(1/N) * (Tanh(t+i*Tau)-1)^(1/N) * Cosh(t+i*Tau)^(2/N)); print "c1:",c1; print "Abs(c1):",Abs(c1);
	//c2 := 1/((Tanh(t - i*Tau)+1)^(1/N) * (Tanh(t-i*Tau)-1)^(1/N) * Cosh(t-i*Tau)^(2/N)); print "c2:",c2; print "Abs(c2):",Abs(c2);
	LHS3 := 1 / (Tanh(t + i*Tau)+1)^(1/N); print "LHS3:",LHS3;
	RHS3 := Cosh(t + i*Tau)^(2/N) * (1 - Tanh(t + i*Tau))^(1/N); print "RHS3:",RHS3;
	print "LHS3/RHS3:",LHS3/RHS3;
	LHS4 := 1 / (Tanh(t - i*Tau)+1)^(1/N); print "LHS4:",LHS4;
	RHS4 := Cosh(t - i*Tau)^(2/N) * (1 - Tanh(t - i*Tau))^(1/N); print "RHS4:",RHS4;
	print "LHS4/RHS4:",LHS4/RHS4;
	"##########";
	LHS5 := 1/((1-Tanh(t))^2 *(1+Tanh(t))^2); print "LHS5:",LHS5;
	RHS5 := Cosh(t)^4; print "RHS5:",RHS5;
	
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



intrinsic RS_SEPMTest( s::RngMPolElt : Prec := 20 ) -> BoolElt, Mtrx
{ Tests if period matrices of superelliptic curves are equivalent }
	print "Computing PM1 with general algorithm:";
	time PM1 := RS_PeriodMatrix(s:Prec:=Prec,LocPar:="x");
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




intrinsic RS_SEAJMTest_1( SEC::SECurve : Ht := 100 ) -> RngIntElt
{ (weak) test for Abel-Jacobi map of superelliptic curves }
	C<i> := SEC`ComplexField;
	R := SEC`RealField;
	x1 := C!Random([-Ht..Ht])/Random([1..Ht]) + i*C!Random([-Ht..Ht])/Random([1..Ht]);
	x2 := C!Random([-Ht..Ht])/Random([1..Ht]) + i*C!Random([-Ht..Ht])/Random([1..Ht]);
	F1 := RS_Fiber(SEC`DefiningPolynomial,x1); print "F1:",F1;
	F2 := RS_Fiber(SEC`DefiningPolynomial,x2); print "F2:",F2;
	Points := < >; Coeffs := [];
	for j in [1..SEC`Degree[1]] do
		Append(~Points,[x1,F1[j]]); Append(~Coeffs,1);
		Append(~Points,[x2,F2[j]]); Append(~Coeffs,-1);
	end for;
	D := RS_SEDivisor(Points,Coeffs,SEC:Check:=true);	
	return SEC`AbelJacobi(D);
end intrinsic;


intrinsic RS_SEAJMTest_2( SEC::SECurve : Ht := 100 ) -> RngIntElt
{ (strong) test for Abel-Jacobi map of superelliptic curves }
	C<i> := SEC`ComplexField; Cxy<x,y> := PolynomialRing(C,2);
	R := SEC`RealField;
	y1 := C!Random([-Ht..Ht])/Random([1..Ht]) + i*C!Random([-Ht..Ht])/Random([1..Ht]); //print "y1:",y1;
	y2 := C!Random([-Ht..Ht])/Random([1..Ht]) + i*C!Random([-Ht..Ht])/Random([1..Ht]); //print "y2:",y2;
	f := Cxy!SEC`DefiningPolynomial;
	g1 := UnivariatePolynomial(Evaluate(f,2,y1)); //print "g1:",g1;
	g2 := UnivariatePolynomial(Evaluate(f,2,y2)); //print "g2:",g2;
	F1 := Roots(g1); //print "F1:",F1;
	F2 := Roots(g2); //print "F2:",F2;
	Points := < >; Coeffs := [];
	for j in [1..SEC`Degree[2]] do
		Append(~Points,[F1[j][1],y1]); Append(~Coeffs,1);
		Append(~Points,[F2[j][1],y2]); Append(~Coeffs,-1);
	end for;
	D := RS_SEDivisor(Points,Coeffs,SEC:Check:=true);
	//print "D:",D;	
	return SEC`AbelJacobi(D);
end intrinsic;


intrinsic RS_SEAJMLTTest_2( N::RngIntElt, d::RngIntElt : Prec:=20, t := 10 ) -> RngIntElt
{ Long time test }
	for j in [1..t] do
		s := RS_RandomSuperellipticCurve(N,d);
		SEC := RS_SECurve(s:Prec:=Prec);
		for k in [1..t] do
			V := RS_SEAJMTest_2(SEC);
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


intrinsic RS_SEAJMTest_Inf( SEC::SECurve )
{ Computes AJM for a principal zero divisor }
	"#####################  TEST ####################";
	delta := SEC`Degree[3];
	assert delta gt 1;
	m := SEC`Degree[1];
	n := SEC`Degree[2];
	delta,a,b := Xgcd(m,n);
	while a ge 0 do
		a -:= n;
		b +:= m;
	end while;
	assert a*m + b*n eq delta;
	M := Round(m/delta);
	N := Round(n/delta);
	C<i> := SEC`ComplexField;
	Ct<t> := PolynomialRing(C);


	// For each points at infinity compute AJM of D = [ r - zeta^k ]
	for k in [0..delta-1] do
		r := Exp(2*SEC`Pi*i*k/delta);
		g_t := &*[ (1 - x*(r+t)^b*(t^M)) : x in SEC`BranchPoints ] - (r+t)^delta;
		g_C := Coefficients(g_t);
		assert Degree(g_t) in [n*(M+b),(n-1)*(M+b)];
		ord := Min([ j : j in [1..Degree(g_t)+1] | Abs(g_C[j]) gt SEC`Error ]) - 1; print "ord:",ord;
		Rts_g_t := RS_Roots(g_t);
		print "Rts:",Rts_g_t;
		//print "Rts_g_t:",Rts_g_t;
		Points := < >; Coeffs := [];
		for j in [1..Degree(g_t)] do
			t_j := Rts_g_t[j];
			if Abs(t_j) gt SEC`Error then
				x_j := 1/(r+t_j)^b/(t_j^M);
				y_j := (r+t_j)^a/((t_j^N)*SEC`LeadingCoeff);
				Append(~Points,[x_j,y_j]);
				Append(~Coeffs,1);
			else
				Append(~Points,[0,k+1,0]);
				Append(~Coeffs,ord);
			end if;
		end for;
		//assert &+Coeffs in [n*b+n*M,n*b+n*M-b-M];
		if Degree(g_t) eq (n-1)*(M+b) then
			"OK";
			Append(~Points,[0,0]);
			Append(~Coeffs,-N*m+M+b);
		else
			Y0 := RS_PrincipalBranch(SEC,Zero(C):Fiber);
			for j in [1..m] do
				Append(~Points,[0,Y0[j]]);
				Append(~Coeffs,-N);
			end for;
		end if;
		for j in [1..n] do
			Append(~Points,[SEC`BranchPoints[j],0]);
			Append(~Coeffs,-b);
		end for;
		D := RS_SEDivisor(Points,Coeffs,SEC:Check:=false);
		//print "D:",D;
		//assert D`Degree eq 0;
		V := SEC`AbelJacobi(D);
		/*print "############### COMPARE #################";
		print "V_test:",V;
		print "-V_real:",SEC`AJM_InftyPoints[k+1];
		print "lol:",V+SEC`AJM_InftyPoints[k+1];*/
		//W :=  b * &+[ V : V in SEC`AJM_RamPoints ]; print "W",W;		
		//V -:= ChangeRing(W,R);
		Z := Matrix(SEC`RealField,2*SEC`Genus,1,[ z - Round(z) : z in Eltseq(V) ]); print "Z:",Z;
		print "Poly:",SEC`DefiningPolynomial;
		assert &and[ Abs(z) lt 10^-10 : z in Eltseq(Z) ];
	end for;
end intrinsic;


intrinsic RS_SEAJMLTTest_Inf( N::RngIntElt, d::RngIntElt : Prec:=20, t := 10 ) -> RngIntElt
{ Long time test }
	for j in [1..t] do
		s := RS_RandomSuperellipticCurve(N,d);
		SEC := RS_SECurve(s:Prec:=Prec);
		for k in [1..t] do
			RS_SEAJMTest_Inf(SEC);
		end for;
	end for;
	return true;
end intrinsic;

intrinsic RS_TestForArb(SEC:SECurve) -> RngIntElt
{ - }
	C<i> := SEC`ComplexField; R := SEC`RealField;
	X := [5, 8 , -9, 8, 1 , 7, -4 , 3 , -5, 0, 0 , -9, 4 , 4, 4, 7 , 10, 4, -3, -9];
	Y := [9, -5, 1 , 0, -4, 8, -10, -8, 0 , 6, -8, -9, -6, 0, 5, -7, 4 , 7, 4 , 5];
	U := [ (X[j]+i*Y[j])/(Sqrt(R!2)) : j in [1..#X] ];

	

/*

  for (j = 0; j < dmax; j++)
        {
            arb_set_si(acb_realref(u + j), x[j]);
            arb_set_si(acb_imagref(u + j), y[j]);
            acb_div_arb(u + j, u + j, s2, prec);
        }
        arb_clear(s2);

        for (j = 0; j < nref; j++)
        {
            if (arb_set_str(acb_realref(ref + j), ref_r[j], prec) ||
                    arb_set_str(acb_imagref(ref + j), ref_i[j], prec) )
            {
                flint_printf("error while setting ref[%ld] <- %s+I*%s\n",
                        j, ref_r[j], ref_i[j]);
                abort();
            }

        }\;
*/

	return 0;
end intrinsic;

	



