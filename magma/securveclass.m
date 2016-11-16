import "superelliptic.m": RS_SEIntegrate, RS_NRootAC, RS_SEMST, RS_MakeCCVectors, RS_SETau3, RS_SEInfTau;

C_20<i> := ComplexField(20);
C_Pi := Real(Pi(C_20));

// Riemann surfaces as type RieSrf defined here 
declare type SECurve;

declare attributes SECurve: DefiningPolynomial, Genus, Degree, BranchPoints, MST,  Tau, HomologyGroup, ComplexField, HolomorphicDifferentials, Prec, Theta, SmallPeriodMatrix, BigPeriodMatrix, Pi, PeriodMatrix, ReductionMatrix, Abscissas, Weights, StepLength, AJWeierstrass, AbelJacobi, TreeMatrix, Error,
ElementaryIntegrals, skd_Matrix, IntersectionMatrix, ABtoC, Basepoint, Zetas, PM_B_Inv;

// Constructor
intrinsic RS_SECurve( f::RngMPolElt : Prec := -1 ) -> SECurve
{ Creates the Riemann surface defined by f(x,y) = 0 and computes its function field and genus }
	
	// Requirements
	N := Degree(f,2); d := Degree(f,1);
	//require &and[N ge 2, d ge 2, N*d ge 6] : "Genus has to be greater than 0.";
	require &and[N ge 2, d ge 3] : "Degrees not supported."; 

	// Create object
	SEC := New(SECurve);	
	
	// Defining polynomial
	SEC`DefiningPolynomial := f;

	// Precision
	if Prec lt 20 then
		SEC`Prec := 25;
	else
		SEC`Prec := Prec+5;
	end if;

	// Error
	SEC`Error := 10^-(SEC`Prec-4);

	// Complex field
	CC<i> := ComplexField(SEC`Prec+10);
	SEC`ComplexField := CC; SEC`Pi := Real(Pi(CC));

	// Degree of cover
	SEC`Degree := [N,d];	

	// Branch points
	f_x := ChangeRing(UnivariatePolynomial(Evaluate(f,[Parent(f).1,0])),CC);
	SEC`BranchPoints :=  RS_Roots(f_x);

	// Holomorphic differentials and genus
	SEC`HolomorphicDifferentials := RS_SuperellipticDifferentials(f);

	// Genus
	SEC`Genus := &+[ #DFF_i : DFF_i in SEC`HolomorphicDifferentials ];

	assert 2*SEC`Genus eq (N-1)*(d-1) - Gcd(N,d) + 1; // Check Riemann-Hurwitz

	// Root of unity powers
	Zeta := RS_ReducePrecision(Exp(2*SEC`Pi*i/N),5); Zeta_ := 1;
	ZetaSqrt := RS_ReducePrecision(Exp(SEC`Pi*i/N),5); ZetaSqrt_ := 1;
	ZetaInv := 1/Zeta; ZetaInv_ := 1;
	ZetaPows := []; ZetaSqrtPows := []; ZetaInvPows := [];
	for j in [1..N-1] do
		ZetaSqrt_ *:= ZetaSqrt;
		Append(~ZetaSqrtPows,ZetaSqrt_);
		ZetaInv_ *:= ZetaInv;
		Append(~ZetaInvPows,ZetaInv_);
		Zeta_ *:= Zeta;
		Append(~ZetaPows,Zeta_);
	end for;
	SEC`Zetas := [ZetaPows,ZetaSqrtPows,ZetaInvPows];

	// Choose branch point with smallest real part as base point
	SEC`Basepoint := 1;

	return SEC;
end intrinsic;

// Printing
intrinsic Print( SEC::SECurve )
{ Print Riemann surface }
	print "Superelliptic curve of genus", SEC`Genus ,"defined by: 0 =",SEC`DefiningPolynomial,"and prescribed precision",SEC`Prec;
	print "";
	print "Computed data:";
	print " Complex field:", SEC`ComplexField;
	print " Branch points:", assigned SEC`BranchPoints;
	print " Roots of unity:", assigned SEC`Zetas;
	print " Maximal spanning tree:", assigned SEC`MST;
	print " Tau for DE-integration:", assigned SEC`Tau;
	print " skd_Matrix:", assigned SEC`skd_Matrix;
	print " Intersection matrix:", assigned SEC`IntersectionMatrix;
	print " Basis of holomorphic differentials:", assigned SEC`HolomorphicDifferentials;
	print " Period matrix (C):", assigned SEC`PeriodMatrix;
	print " Big period matrix (A B):",assigned SEC`BigPeriodMatrix;
	print " Small period matrix (B^-1 A):",assigned SEC`SmallPeriodMatrix;
	print " Theta function:",assigned SEC`Theta;
	print " Abel-Jacobi map:", assigned SEC`AbelJacobi;
	print " AJWeierstrass:", assigned SEC`AJWeierstrass;
	print " Basepoint:", assigned SEC`Basepoint;
	print " Elementary integrals:", assigned SEC`ElementaryIntegrals;
	print " Tree matrix:", assigned SEC`TreeMatrix;
	print " Reduction matrix:", assigned SEC`ReductionMatrix;
end intrinsic;


intrinsic RS_ChooseAJPath( SEC::SECurve, P::FldComElt ) -> RngIntElt
{ Choose path such that tau for integration is maximal }
	BestTau := 0; BestInd := 0;
	Tau := SEC`Tau; d := SEC`Degree[2]; Pts := SEC`BranchPoints;
	for j in [1..d] do
		Tau_k := 4.0; // > Pi/2
		for k in [1..d] do
			if k ne j then
				Tau_k := Min(Tau_k,RS_SETau3(P,Pts[j],Pts[k]));
			end if;
		end for;
		if Tau_k gt BestTau then
			BestInd := j;
			BestTau := Tau_k;
		end if;
	end for;

	if BestTau lt SEC`Tau then
		// Update integration parameters
		print "Old tau:",SEC`Tau;
		SEC`Tau := BestTau;
		print "New tau:",SEC`Tau;
		SEC`Abscissas, SEC`Weights, SEC`StepLength := RS_SEIntegrationParameters(Pts,SEC`MST,BestTau,SEC`Degree[1],SEC`ComplexField);
	end if;

	return BestInd;
	
end intrinsic;


intrinsic RS_LatticeReduction( SEC::SECurve, Z::[FldComElt] ) -> RngIntElt
{ Reduce z \in \C^g modulo period matrix (A B) // Really? }
	R := RealField(SEC`Prec); g := SEC`Genus; 
	if not assigned SEC`ReductionMatrix then
		if not assigned SEC`BigPeriodMatrix then
			RS_SEPM(SEC);
		end if;
		PM_AB := SEC`BigPeriodMatrix;
		M := ZeroMatrix(R,2*g,2*g);
		for j in [1..g] do
			for k in [1..g] do
				M[j,k] := Re(PM_AB[j,k]);
				M[j+g,k] := Im(PM_AB[j,k]);
				M[j,k+g] := Re(PM_AB[j,k+g]);	
				M[j+g,k+g] := Im(PM_AB[j,k+g]);
			end for;
		end for;
		SEC`ReductionMatrix := M^(-1);
	end if;

	RProj := SEC`ReductionMatrix;
	Rz := Matrix(R,2*g,1,[ Re(z) : z in Z ] cat [ Im(z) : z in Z ]);
	centerZ := RProj * Rz;
	return centerZ;
end intrinsic;
intrinsic RS_ModPeriodLattice( SEC::SECurve, V::AlgMatElt[FldCom] ) -> SeqEnum[FldComElt]
{ Reduce the vector V \in \C^g / <1,PM>, see van Wamelen's code }
	RS_SEPM(SEC);
	g := SEC`Genus;
	PM := SEC`SmallPeriodMatrix;
	C := BaseRing(PM);
	I_PM := Matrix(g,g,[Im(a) : a in ElementToSequence(PM)]);
	dum := I_PM^(-1)*Matrix(g,1,[Im(zi) : zi in ElementToSequence(V)]);
	v1 := Matrix(C,g,1,[Round(di) : di in ElementToSequence(dum)]);
	V := V - PM*v1;
	return V - Matrix(C,g,1,[Round(Re(di)) : di in ElementToSequence(V)]);
end intrinsic;
intrinsic RS_ModPeriodLattice( SEC::SECurve, V::ModMatFldElt[FldCom] ) -> SeqEnum[FldComElt]
{ Reduce the vector V \in \C^g / <1,PM>, see van Wamelen's code }
	RS_SEPM(SEC);
	g := SEC`Genus;
	PM := SEC`SmallPeriodMatrix;
	C := BaseRing(PM);
	I_PM := Matrix(g,g,[Im(a) : a in ElementToSequence(PM)]);
	dum := I_PM^(-1)*Matrix(g,1,[Im(zi) : zi in ElementToSequence(V)]);
	v1 := Matrix(C,g,1,[Round(di) : di in ElementToSequence(dum)]);
	V := V - PM*v1;
	return V - Matrix(C,g,1,[Round(Re(di)) : di in ElementToSequence(V)]);
end intrinsic;

intrinsic RS_UpperSheet( SEC::SECurve, z::FldComElt : Global := true ) -> FldComElt
{ - }
	Pts := SEC`BranchPoints;
	if Global then
		return <z,(&*[ z-Pts[k] : k in [1..SEC`Degree[2]] ])^(1/SEC`Degree[1])>;
	else
		return <z,(&*[ (z-Pts[k])^(1/SEC`Degree[1]) : k in [1..SEC`Degree[2]] ])>;
	end if;
end intrinsic;
					
intrinsic RS_AJIntegrate( SEC::SECurve, P::Tup, Ind::RngIntElt ) -> SeqEnum[FldComElt]
{ Integrate from P to P_i }
	I := RS_AJIntegrate( SEC, Ind, P );
	return [ -I[j] : j in [1..SEC`Genus] ];
end intrinsic;
intrinsic RS_AJIntegrate( SEC::SECurve, Ind::RngIntElt, P::Tup ) -> SeqEnum[FldComElt]
{ Integrate from P_i to P }
	CC<i> := SEC`ComplexField; Pi := SEC`Pi; CC_0 := Zero(CC); CC_1 := One(CC); g := SEC`Genus;
	Zeta := SEC`Zetas[1][1]; d := SEC`Degree[2]; N := SEC`Degree[1]; Points := SEC`BranchPoints; DFF := &cat[ DFF_i : DFF_i in SEC`HolomorphicDifferentials ];

	a := Points[Ind];
	P_x := P[1];
	P_y := P[2];
	
	// Integrate from left to right
	// Not useful here?
	//assert Re(a) le Re(p);
	
	// Make vector of centers of circumcircles
	CCV := RS_MakeCCVectors(Ind,P_x,Points);
	assert #CCV eq d-1;

	// Factors due to change of variable
	Fact1 := (P_x-a)/2;
	Fact2 := (P_x+a)/(P_x-a);
	
	vprint SE,2 : "################### Next Integral ###################";
	Integral := [ CC_0 : j in [1..g] ];
	
	// Initiate on x = 0, dx = 1
	z := RS_NRootAC(0,CCV,Zeta,N);
	for j in [1..g] do
		w := DFF[j];
		Denom := z^w[2];
		Integral[j] +:= (Fact2^w[1]/Denom);
	end for;

	for t in [1..#SEC`Abscissas] do
		x := SEC`Abscissas[t];

		// Analytic continuation
		z1 := RS_NRootAC(x,CCV,Zeta,N);
		z2 := RS_NRootAC(-x,CCV,Zeta,N);
		
		// End point not a branch point
		c1 := (1-x)^(1/N);
		c2 := (1+x)^(1/N);

		//print "c1:",c1; print "c2:",c2;

		Enum1 := (x + Fact2);
		Enum2 := (-x + Fact2);

		for j in [1..g] do
			w := DFF[j];
			Denom1 := z1^w[2];
			Denom2 := z2^w[2];
			dx := SEC`Weights[t][1] * SEC`Weights[t][2]^(2-(2*w[2]/N));
			Integral[j] +:= ( Enum1^w[1] * c1^w[2] / Denom1 + Enum2^w[1] * c2^w[2] / Denom2) * dx;
		end for;
	end for;

	// Correct sheet?
	
	P_y_AC := ((P_x-a)/2)^(d/N) * RS_NRootAC(1,CCV,Zeta,N); // * Cosh( w(P_x) )^(2/N) * 2^(1/N)

	//print "AJIntegrate0:",Integral;
	//print "P_y:",P_y;
	//print "P_y_AC:",P_y_AC;
	ValArg := Arg(P_y_AC);
	print "Arg(P_y_AC):",ValArg;
	print "Arg P_y:",Arg(P_y);
	k := (N/(2*SEC`Pi)) * ( Arg(P_y) - ValArg );
	print "k:",k;
	k_ := Round(k);
	assert Abs(k - k_) lt 10^-10; // Check: k \in \Z ?
	k := k_ mod N;
	for j in [1..g] do
		w := DFF[j];
		Factor := Fact1^((w[1]+1)-(d*w[2]/N));
		//Integral[j] *:= SEC`Zetas[1][w[2]]^k * SEC`Zetas[2][w[2]] * SEC`StepLength * Factor;
		Integral[j] *:= SEC`Zetas[1][w[2]]^k * SEC`StepLength * Factor;
	end for;
	print "AJIntegrate:",Integral;
	return Matrix(CC,g,1,Integral);
end intrinsic;


function RS_NRootACInfty(x,p,Pts,Zeta,N:long:=false)
// Analytic continuation of y = f(x)^(1/N)  for x -> \infty on one sheet, avoiding branch cuts
	long := true;
	if #Pts eq 1 or long then
		prod := 1;
		for k in [1..#Pts] do
			prod *:= ( 2*p + (x-1) * Pts[k]  )^(1/N);
		end for;
		return prod;
	end if;
end function;

intrinsic RS_AJInfinity( SEC::SECurve ) -> SeqEnum[FldComElt]
{ Integrate from \infty to P_0 }

	RS_ComputeAll(SEC);

	"################## Integrate \infty to P_0 ####################";

	CC<i> := SEC`ComplexField;
	CC_0 := Zero(CC);
	Pts := SEC`BranchPoints;
	Zeta := SEC`Zetas[1][1];	
	N := SEC`Degree[1];
	d := SEC`Degree[2];
	g := SEC`Genus;
	print "Pts:",Pts;

	x_0 := Pts[SEC`Basepoint];

	
	print "x_0:",x_0;

 	MaxIm := Max( [Abs(Im(Pts[k])) : k in [1..d] ]);

	print "MaxIm:",MaxIm;

	p := 2*(Min(Re(x_0)-1,-1) + i * (MaxIm + 1));
	//p := CC!(1.6666666666666666667 - 1.6665612762750131730*i);
	//p := Re(x_0) + 1000*i*(MaxIm+1);
	print "p:",p;

	InfTau := RS_SEInfTau( p, Pts : Lambda := C_Pi/2 );
	print "InfTau:",InfTau;
	assert InfTau ge SEC`Tau;
	// TODO Better Tau?

	P := RS_UpperSheet(SEC,p);

	print "P:",P;

	// Integrate from P_0 to P
	V1 := RS_AJIntegrate(SEC,1,P);

	// Integrate from P to \infty
	Fact := 2*p;
	print "Fact:",Fact;
	Integral := [ CC_0 : j in [1..g] ];
	DFF := &cat[ DFF_i : DFF_i in SEC`HolomorphicDifferentials ];

	// Initiate on x = 0, dx = 1
	z := RS_NRootACInfty(0,p,Pts,Zeta,N);
	for j in [1..g] do
		w := DFF[j];
		Denom := z^w[2];
		Integral[j] +:= (1/Denom);
	end for;

	for t in [1..#SEC`Abscissas] do
		x := SEC`Abscissas[t];

		// Analytic continuation
		z1 := RS_NRootACInfty(x,p,Pts,Zeta,N);
		z2 := RS_NRootACInfty(-x,p,Pts,Zeta,N);
		print "z1:",z1;
		Enum1 := (1 - x);
		Enum2 := (1 + x);

		for j in [1..g] do
			w := DFF[j];
			pow := d*w[2]/N - w[1] - 1;
			assert pow ge 0;
			Denom1 := z1^w[2];
			Denom2 := z2^w[2];
			dx := SEC`Weights[t][1];
			//dx := SEC`Weights[t][1] * SEC`Weights[t][2]^2;
			//Integral[j] +:= ( Enum1^w[1] * c1^w[2] / Denom1 + Enum2^w[1] * c2^w[2] / Denom2) * dx;
			Integral[j] +:= ( (Enum1^pow * Enum2 / Denom1) + (Enum2^pow * Enum1 / Denom2) ) * dx;
		end for;
	end for;


	P_y_AC := RS_NRootACInfty(-1,p,Pts,Zeta,N); // * 2^(d/N)
	//print "P_y_AC:",P_y_AC;
	print "Arg(P_y_AC):",Arg(P_y_AC);
	//print "P_y:",P[2];
	print "Arg(P_y):",Arg(P[2]);	

	k := (N/(2*SEC`Pi)) * ( Arg(P[2]) - Arg(P_y_AC) );
	print "k:",k;
	k_ := Round(k);
	assert Abs(k - k_) lt 10^-10; // Check: k \in \Z ?
	k := k_ mod N;
	
	for j in [1..g] do
		w := DFF[j];
		Integral[j] *:= SEC`Zetas[1][w[2]]^k * Fact^(w[1]+1);
	end for;

	print "Integral_x0toP",V1;
	print "Integral_PtoInf:",Integral;
	return Integral;
	// I_{\infty} = (-1)(I_{P_0,P} + I_{P,\infty})
	Integral_Inftox_0 := -(Matrix(BaseRing(SEC`SmallPeriodMatrix),g,1,Integral) + V1);

	print "Integral_Inftox0:",Integral_Inftox_0;

	V := RS_ModPeriodLattice(SEC,ColumnSubmatrix(Integral_Inftox_0,1,1));

	return V;
end intrinsic;


intrinsic RS_SymplecticInverse( S::Mtrx ) -> Mtrx
{ Computes the inverse of a symplectic matrix }
	1 eq 1;
end intrinsic;

intrinsic RS_AbelJacobi( SEC::SECurve : Recompute := false )
{ Compute the Abel-Jacobi map of the point P = (x,y) }

	if Recompute or not assigned SEC`AbelJacobi then

		// Complex Field
		C<i> := SEC`ComplexField;

		// Scaling periods to half-periods // False?	
		DFF := SEC`HolomorphicDifferentials; SV := [];
		for j in [1..#DFF] do
			for k in [1..#DFF[j]] do
				//Append(~SV,SEC`Zetas[2][DFF[j][k][2]]/(SEC`Zetas[2][DFF[j][k][2]]-SEC`Zetas[3][DFF[j][k][2]]));
				Append(~SV,1/(1-(SEC`Zetas[3][DFF[j][k][2]])));
			end for;
		end for;
		SV := RS_Vector(SV);
	
		// Compute period matrix
		RS_SEPM(SEC:Recompute:=Recompute);
		ABtoC_Inv := SEC`ABtoC^(-1);
		// Compute 'map' of the tree
		RS_TreeMatrix(SEC:Recompute:=Recompute);
		// 'A-periods'
		SEC`PM_B_Inv := ColumnSubmatrix(SEC`BigPeriodMatrix,1,SEC`Genus)^(-1);
		// AJ between Weierstrass point
		RS_AJWeierstrass(SEC:Recompute:=Recompute);

		

		// Define Abel-Jacobi map
		AbelJacobi := function ( P )
		// Computes the vector of integrals of holomorphic differentials along a path from P_0 to P
			"########### Start AJ #############";
			

			if #P eq 1 then
				x_P := C!P[1];
				y_P := RS_UpperSheet(SEC,C!x_P)[2];
			elif #P eq 2 then
				x_P := C!P[1];
				y_P := C!P[2];
			else
				error Error(".");
			end if;

			// Check if P is a branch point
			Dist, Ind := RS_Distance(x_P,SEC`BranchPoints);

			print "Dist:",Dist; print "Ind:",Ind;	
			//V := RS_ZeroVector(C,SEC`Genus);
			if Dist lt SEC`Error then
				/*
				TreePath := SEC`TreeMatrix[Ind];
				print "TreePath:",TreePath;
				print "SV:",SV;
				W := RS_AJWeierstrass(SEC: Ind:=Ind);
				print "W:",W;
				for j in [1..SEC`Degree[2]-1] do
					//V +:= TreePath[j] * RS_Vector(SEC`ElementaryIntegrals[j]);
					V +:= TreePath[j] * SV;
				end for;
				print "V:",V;
				//V := PM_B_Inv * V;
				return V;
				return RS_ModPeriodLattice(SEC,V);
				print "V:",V;
				print "N*V:",SEC`Degree[1]*V;
				*/
				//return ABtoC_Inv * V;
				
				return ColumnSubmatrix(SEC`AJWeierstrass,Ind,1);
			end if;

			Ind := RS_ChooseAJPath( SEC, x_P );
			

			// Integrate from P_i to P
			I := RS_AJIntegrate(SEC,Ind,<x_P,y_P>);

			// AJM of Weierstrass points
			W := ColumnSubmatrix(SEC`AJWeierstrass,Ind,1);

			// Sum
			V := I+W;

			// Reduce mod lattice
			return RS_ModPeriodLattice(SEC,V);

			// return RS_LatticeReduction(SEC,V);
			
		end function;

		SEC`AbelJacobi := AbelJacobi;
	end if;	

end intrinsic;


intrinsic RS_TreeMatrix( SEC::SECurve : Recompute := false  )
{ Compute a matrix with paths from x_0 -> x_i for all branch points }

	if not assigned SEC`TreeMatrix or Recompute then

	// Compute period matrix
	RS_SEPM(SEC);
	
  	Pts := SEC`BranchPoints;
 	g := SEC`Genus; N := SEC`Degree[1]; d := SEC`Degree[2];
 
	TM := ZeroMatrix(Integers(),d,d-1);
	
	Taken := [ 0 : j in [1..d] ];
	Tree := SEC`MST;
	P_0 := Tree[1][1];
	Taken[P_0] := 1;
	
	for j in [1..#Tree] do
		if Taken[Tree[j][1]] eq 1 then
			PStart := Tree[j][1];
			PEnd := Tree[j][2];
			Taken[Tree[j][2]] := 1;
		else
			PStart := Tree[j][2];
			PEnd := Tree[j][1];
			Taken[Tree[j][1]] := 1;
		end if;
		TM[PEnd] := TM[PStart];

      		TM[PEnd,j] := 1;
	end for;

	P_0 := SEC`Basepoint; // Shift by real basepoint
	P_0P_0 := TM[P_0];
	for j in [1..d] do
		TM[j] -:= P_0P_0;
	end for;

	SEC`TreeMatrix := TM;

	end if;

end intrinsic;


intrinsic RS_AJWeierstrass( SEC::SECurve : Recompute := false )
{ Compute image of branch points under Abel-Jacobi map }

	if not assigned SEC`AJWeierstrass or Recompute then

	// Compute period matrix
	RS_SEPM(SEC); C<i> := CoefficientRing(SEC`BigPeriodMatrix);
	
  	Pts := SEC`BranchPoints;
 	g := SEC`Genus; N := SEC`Degree[1]; d := SEC`Degree[2];
 
	//P0Pi := ZeroMatrix(Rationals(),d,(N-1)*(d-1));
	P0Pi := Matrix(C,d,(N-1)*(d-1),[]);
	
	Taken := [ 0 : j in [1..d] ];
	Tree := SEC`MST;
	P0 := Tree[1][1];
	Taken[P0] := 1;
	assert P0 eq SEC`Basepoint;
	for j in [1..#Tree] do
		if Taken[Tree[j][1]] eq 1 then
			PStart := Tree[j][1];
			PEnd := Tree[j][2];
			Taken[Tree[j][2]] := 1;
		else
			PStart := Tree[j][2];
			PEnd := Tree[j][1];
			Taken[Tree[j][1]] := 1;
		end if;
		P0Pi[PEnd] := P0Pi[PStart];
      		P0Pi[PEnd,(N-1)*(j-1)+1] := 1/N;
	end for;

	P0 := SEC`Basepoint; // Shift by real basepoint
	P0P0 := P0Pi[P0];
	for j in [1..d] do
		P0Pi[j] -:= P0P0;
	end for;

	// g x d Matrix with image of AJM between Weierstrass points
	M1 := ChangeRing(SEC`ABtoC^(-1),C) * Transpose(P0Pi); // = (d-1)(N-1)x(d) = (d-1)(N-1)x(d-1)(N-1) x (d-1)(N-1)x(d)
	M2 := Submatrix(M1,[1..2*g],[1..d]);
	V := HorizontalJoin(DiagonalMatrix(C,[ One(C) : j in [1..g]]),SEC`SmallPeriodMatrix) * M2; // (g)x(d) = (g)x(2g) x (2g) x (d)
	//V := SEC`PM_B_Inv * SEC`BigPeriodMatrix * ChangeRing(SEC`ABtoC^(-1),C) * Transpose(P0Pi);
	
	assert Nrows(V) eq SEC`Genus;

	SEC`AJWeierstrass := V;
	
	end if;
end intrinsic;



intrinsic RS_SEPM( SEC::SECurve : Recompute := false, Small := true )
{ Computes period matrices associated to the superelliptic curve defined by f to precision of SEC }

	if Recompute or &and[Small,not assigned SEC`SmallPeriodMatrix] or &and[not Small, not assigned SEC`BigPeriodMatrix] then 

	// Branch points
	Points := SEC`BranchPoints;
	
	// Degrees
	N := SEC`Degree[1]; d := SEC`Degree[2]; Err := 10^-(SEC`Prec+1);

	// Complex fields and constants
	C_<i> := ComplexField(SEC`Prec); Zeta := SEC`Zetas[1][1]; vprint SE,1 : "Zeta:",C_20!Zeta;
	CC<i> := SEC`ComplexField; CC_0 := Zero(CC);

	
	// Maximal spanning tree w.r.t holomorphicy
	vprint SE,1 : "Constructing maximal spanning tree...";
	t := Cputime();
	SEC`MST, SEC`Tau := RS_SEMST(Points);
	Cputime(t);
	vprint SE,2 : "MST_Edges:",SEC`MST;
	vprint SE,2 : "#MST_Edges:",#SEC`MST;
	vprint SE,2 : "MST_Tau:",SEC`Tau;

	// Holomorphic differentials
	DFF := SEC`HolomorphicDifferentials;

	// Genus
	g := SEC`Genus;
	vprint SE,1 : "Genus:",g;

	// Integration parameters
	vprint SE,1 : "Computing Integration parameters...";
	t := Cputime();
	SEC`Abscissas, SEC`Weights, SEC`StepLength := RS_SEIntegrationParameters(Points,SEC`MST,SEC`Tau,N,SEC`ComplexField);
	vprint SE,2 : "#Abscissas:",#SEC`Abscissas;
	vprint SE,3 : "Abscissas:",SEC`Abscissas;
	vprint SE,2 : "#Weights:",#SEC`Weights;
	vprint SE,3 : "Weights:",SEC`Weights;
	vprint SE,2 : "StepLength:",SEC`StepLength;
	vprint SEP,1 : "Precision Abscissas:",Precision(SEC`Abscissas[1]);
	vprint SEP,1 : "Precision Weights:",Precision(SEC`Weights[1][1]);
	vprint SEP,1 : "Precision StepLength:",Precision(SEC`StepLength);
	Cputime(t);

	// Integrals
	vprint SE,1 : "Integrating...";
	t := Cputime();
	Integrals := [];
	SEC`ElementaryIntegrals := [];
	for k in [1..d-1] do
		vprint SE,2 : "Integrating edge nr.",k;
		I, EI := RS_SEIntegrate( SEC`MST[k],Points,DFF,SEC`Abscissas,SEC`Weights,SEC`StepLength,d,N,SEC`Zetas );
		Append(~SEC`ElementaryIntegrals,EI);
		Append(~Integrals,I);
	end for;
	Cputime(t);


	// Matrix of periods
	PM_C := ZeroMatrix(CC,g,(d-1)*(N-1));
	for k in [1..d-1] do
		for l in [1..N-1] do
			Ind2 := (N-1)*(k-1) + l;
			for j in [1..g] do
				PM_C[j][Ind2] := Integrals[k][j][l];
			end for;
		end for;
	end for;
	vprint SE,3: "Integrals:",PM_C;
	SEC`PeriodMatrix := PM_C;

	// Intersection matrix
	vprint SE,1 : "Computing skd-Matrix...";
	t := Cputime();
	skd_Matrix := [ [] : j in [1..d-1] ];
	for j in [1..d-1] do
		skd_Matrix[j][j] := <0,0,0>;
		for l in [j+1..d-1] do
			skd := RS_SEIntersection(SEC`MST[j],SEC`MST[l],Points,N,Zeta);
			vprint SE,2: "<s,k,d>-triple of edges",SEC`MST[j],"and",SEC`MST[l],":",skd;
			skd_Matrix[j][l] := <skd[1],skd[2] mod N,skd[3]>;
			skd_Matrix[l][j] := <-skd[1],skd[2] mod N,skd[3]>;
		end for;
	end for;
	Cputime(t);
	SEC`skd_Matrix := skd_Matrix;
	vprint SE,2: "skd_matrix:",skd_Matrix;

	// Building block matrices
	Blocks := [];
	for j in [1..d-1] do
		for l in [1..d-1] do
			Block := ZeroMatrix(Integers(),N-1,N-1);
			//Block := ZeroMatrix(Integers(),N,N);
			if j eq l then
				// Intersection with self-lifts
				if N gt 2 then
					for Ind1 in [1..N-1] do
					//for Ind1 in [1..N] do
						for Ind2 in [Ind1+1..N-1] do
						//for Ind2 in [Ind1+1..N] do
							if Ind1 + 1 eq Ind2 then
								Block[Ind1][Ind2] := 1;
								Block[Ind2][Ind1] := -1;
							end if;
						end for;
					end for;
					//Block[1][N] := -1;
					//Block[N][1] := 1;
				end if;
			else
				s := skd_Matrix[j][l][1];
				if s ne 0 then
					k := skd_Matrix[j][l][2];
					dir := skd_Matrix[j][l][3];
					for Ind1 in [1..N-1] do
					//for Ind1 in [1..N] do
						for Ind2 in [1..N-1] do
						//for Ind2 in [1..N] do
							Shift := (Ind2 - Ind1 - k) mod N;
							if Shift eq 0 then
								Block[Ind1][Ind2] := s;
							elif Shift eq dir then
								Block[Ind1][Ind2] := -s;					
							end if;
						end for;
					end for;
				end if;
			end if;
			if j gt l then
				Append(~Blocks,Transpose(Block));
			else
				Append(~Blocks,Block);
			end if;
		end for;
	end for;
	//print "Blocks:",Blocks;
	SEC`IntersectionMatrix := BlockMatrix(d-1,d-1,Blocks);
	vprint SE,3: "(Real) Intersection matrix:",SEC`IntersectionMatrix;
	vprint SE,2: "Rank(IntersectionMatrix):",Rank(SEC`IntersectionMatrix);
	assert Rank(SEC`IntersectionMatrix) eq 2*g;	
	
	// Symplectic reduction of intersection matrix
	vprint SE,1 : "Performing symplectic reduction...";
	t := Cputime();
	CF, ABtoC := RS_SymplecticBasis(SEC`IntersectionMatrix);
	SEC`ABtoC := Transpose(ABtoC);
	ABtoC_CC := ChangeRing(SEC`ABtoC,CC);
	Cputime(t);
	vprint SE,3: "ST:",SEC`ABtoC;
	vprint SE,3: "CF:",CF;
	
	
	
	vprint SE,1 : "Matrix multiplication 1...";
	t := Cputime();
	PMAPMB := PM_C * ABtoC_CC;
	Cputime(t);
	
	PM_A := ColumnSubmatrixRange(PMAPMB,1,g);
	PM_B := ColumnSubmatrixRange(PMAPMB,g+1,2*g);

	vprint SE,3 : "Dependent columns:",ColumnSubmatrixRange(PMAPMB,2*g+1,Nrows(ABtoC_CC));
	vprint SE,3 : "PM_A:",PM_A;
	vprint SE,3 : "PM_B:",PM_B;
	//vprint SE,2 : "Det(PM_A):",Determinant(PM_A);
	//vprint SE,2 : "Det(PM_B):",Determinant(PM_B);
	//vprint SE,3 : "PM_B^-1:",PM_B^-1;
	

	// Compute big period matrix
	PM := HorizontalJoin(PM_B,PM_A) + ZeroMatrix(C_,g,2*g); // Dimensions???
	vprint SE,1 : "Period matrix:";
	SEC`BigPeriodMatrix := PM;


	if Small then
		vprint SE,1 : "Matrix inversion...";
		t := Cputime();
		PM_BInv := PM_B^-1;
		Cputime(t);
		vprint SE,1 : "Matrix multiplication 2...";
		t := Cputime();
		PM := PM_BInv * PM_A;
		Cputime(t);
	
		// Smoothing out entries
		vprint SE,1 : "Smoothing out entries...";
		t := Cputime();
		for j in [1..g] do
			for k in [1..g] do
				if Abs(Re(PM[j][k])) lt Err then
					PM[j][k] := Zero(C_) + Im(PM[j][k])*i;
				end if;
				if Abs(Im(PM[j][k])) lt Err then
					PM[j][k] := Re(PM[j][k]) + Zero(C_)*i;		
				end if;
			end for;
		end for;
		Cputime(t);
		PM := ChangeRing(PM,C_);
		vprint SE,3 : "PM:",PM;


		// Testing for symmetry of the period matrix
		vprint SE,1 : "Testing symmetry...";
		t := Cputime();
		MaxSymDiff := 0;
		for j in [1..g] do
			for k in [j+1..g] do
				MaxSymDiff := Max(MaxSymDiff,Abs(PM[j][k] - PM[k][j]));
			end for;
		end for;
		Cputime(t);
		vprint SE,1 : "Maximal symmetry deviation:",MaxSymDiff;
	
		if MaxSymDiff ge 10^-10 then
			error Error("Period matrix not symmetric: Computation failed.");
		end if;	
		
		// Testing positive definiteness of the imaginary part of the period matrix
		vprint SE,1 : "Testing positive definite...";
		t := Cputime();
		PM_Im := ZeroMatrix(RealField(SEC`Prec),g,g);
		for j in [1..g] do
			for k in [j..g] do
				PM_Im[j][k] := Real(Im(PM[j][k]));
				PM_Im[k][j] := PM_Im[j][k];
			end for;
		end for;
		assert IsPositiveDefinite(PM_Im);
		Cputime(t);
		SEC`SmallPeriodMatrix := PM;
	end if;

	end if;
end intrinsic;


// Some functionality
intrinsic RS_DeleteAll( SEC::SECurve )
{ Deletes all computed data }

	delete SEC`MST;
	delete SEC`Tau;
	delete SEC`HomologyGroup;
	delete SEC`SmallPeriodMatrix;
	delete SEC`BigPeriodMatrix;
	delete SEC`AbelJacobi;
	delete SEC`Theta;
	delete SEC`skd_Matrix;
	delete SEC`IntersectionMatrix;
	delete SEC`ABtoC;
	delete SEC`ElementaryIntegrals;

end intrinsic;

intrinsic RS_ComputeAll( SEC::SECurve : Recompute := false )
{ Computes all data }

	if Recompute then
		RS_DeleteAll(SEC);
	end if;

	RS_AbelJacobi(SEC:Recompute:=Recompute);

end intrinsic;


