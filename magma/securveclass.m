// Import global settings
import "superelliptic.m": RS_SEIntegrate, RS_NRootAC, RS_SEMST;
import "comparefunctions.m": RS_CompareByFirstComplexEntry,RS_CompareFldComElt;
import "globalprecision.m": RS_Config;
import "pathmethods.m": RS_ChoiceOfPathStar;



// Riemann surfaces as type RieSrf defined here 
declare type SECurve;

declare attributes SECurve: DefiningPolynomial, Genus, Degree, BranchPoints, MST,  MST_Tau, HomologyGroup, HolomorphicDifferentials, Prec, Theta, SmallPeriodMatrix, BigPeriodMatrix, AbelJacobi,
			    skd_Matrix, IntersectionMatrix, SymplecticTransformation;

// Constructor
intrinsic RS_SECurve( f::RngMPolElt : Prec := -1 ) -> SECurve
{ Creates the Riemann surface defined by f(x,y) = 0 and computes its function field and genus }
	
	// Requirements
	N := Degree(f,2); d := Degree(f,1);
	//require &and[N ge 2, d ge 2, N*d ge 6] : "Genus has to be greater than 0.";
	require &and[N ge 2, d ge 3] : "Degrees not supported."; 
	require Gcd(f,Derivative(f,1)) eq 1 : "Curve has to be non-singular.";
	require BaseRing(Parent(f)) eq Rationals() : "Polynomial has to be defined over \Q.";

	// Create object
	SEC := New(SECurve);	
	
	// Defining polynomial
	SEC`DefiningPolynomial := f;

	// Degree of cover
	SEC`Degree := N;

	// Precision
	if Prec lt 20 then
		SEC`Prec := RS_GetGlobalPrecision();
	else
		SEC`Prec := Prec;
	end if;

	// Branch points
	SEC`BranchPoints := Sort(RS_DiscriminantPoints(f,Parent(f).2:Prec:=SEC`Prec),RS_CompareFldComElt);

	// Holomorphic differentials and genus
	SEC`HolomorphicDifferentials := RS_SuperellipticDifferentials(f);

	// Genus
	SEC`Genus := #SEC`HolomorphicDifferentials;

	assert 2*SEC`Genus eq (N-1)*(d-1) - Gcd(N,d) + 1; // Check Riemann-Hurwitz

	return SEC;
end intrinsic;

// Printing
intrinsic Print( SEC::SECurve )
{ Print Riemann surface }
	print "Superelliptic curve of genus", SEC`Genus ,"defined by: 0 =",SEC`DefiningPolynomial,"and prescribed precision",SEC`Prec;
	print "";
	print "Computed data:";
	print " Branch points:", assigned SEC`BranchPoints;
	print " Maximal spanning tree:", assigned SEC`MST;
	print " Tau for DE-integration:", assigned SEC`MST_Tau;
	print " skd_Matrix:", assigned SEC`skd_Matrix;
	print " Intersection matrix:", assigned SEC`IntersectionMatrix;
	print " Basis of holomorphic differentials:", assigned SEC`HolomorphicDifferentials;
	print " Big period matrix:",assigned SEC`BigPeriodMatrix;
	print " Small period matrix:",assigned SEC`SmallPeriodMatrix;
	print " Theta function:",assigned SEC`Theta;
	print " Abel-Jacobi map:", assigned SEC`AbelJacobi;
end intrinsic;

/*
intrinsic RS_DeleteAll( X::RieSrf )
{ Deletes all computed data }

	delete X`DiscriminantPoints;
	delete X`ExactDiscriminantPoints;
	delete X`EmbeddingMap;
	delete X`Basepoint;
	//delete X`PathMethod;
	delete X`PathPieces;
	delete X`IndexPathLists;
	delete X`FundamentalGroup;
	delete X`BranchPoints;
	delete X`LocalMonodromy;
	delete X`MonodromyGroup;
	delete X`HomologyGroup;
	delete X`BasisOfDifferentialsFirstKind;
	//delete X`IntMethod;
	delete X`PeriodMatrix;
	delete X`Theta;
	delete X`AbelJacobi;
	delete X`PuiseuxSeries;
	

end intrinsic;

intrinsic RS_ComputeAll( X::RieSrf : Recompute := false )
{ Computes all data }

	if Recompute then
		RS_DeleteAll(X);
	end if;

	RS_DiscriminantPoints(X : Exact);
	RS_PathPieces(X);
	RS_Monodromy(X);
	RS_FundamentalGroup(X);
	RS_HomologyGroup(X);
	RS_BasisOfDifferentialsFirstKind(X);
	RS_PeriodMatrix(X: IntMethod := X`IntMethod);
	RS_Theta(X);
	RS_AbelJacobi(X);
	RS_PuiseuxSeries(X);

end intrinsic;
*/


intrinsic RS_SEPM( SEC::SECurve : Small := true )
{ Computes period matrices associated to the superelliptic curve defined by f to precision of SEC }

	// Branch points
	Points := SEC`BranchPoints;

	// Degrees
	N := SEC`Degree; d := #Points;

	// Complex fields and constants
	C_<i> := ComplexField(SEC`Prec); C_20<i> := RS_GetGlobalComplexField_20(); RS_SetGlobalZeta(N); Zeta := RS_GetGlobalZeta(); vprint SE,1 : "Zeta:",C_20!Zeta;
	CC<i> := RS_GetGlobalComplexField_Max(); CC_0 := Zero(CC); Pi := RS_GetGlobalPi(); 

	// Speial case for hyperelliptic curves of even degree
	if N eq 2 and d mod 2 eq 0 then
		Prune(~Points);
		d -:= 1;
	end if;
	
	// Maximal spanning tree w.r.t holomorphicy
	vprint SE,1 : "Constructing maximal spanning tree...";
	t := Cputime();
	SEC`MST, SEC`MST_Tau := RS_SEMST(Points);
	Cputime(t);
	vprint SE,3 : "MST_Edges:",SEC`MST_Edges;
	vprint SE,2 : "#MST_Edges:",#SEC`MST_Edges;
	vprint SE,2 : "MST_Tau:",SEC`MST_Tau;

	// Holomorphic differentials
	DFF := SEC`HolomorphicDifferentials;

	// Genus
	g := SEC`Genus;
	vprint SE,1 : "Genus:",g;

	// Integration parameters
	vprint SE,1 : "Computing Integration parameters...";
	t := Cputime();
	Abscissas, Weights, StepLength := RS_SEIntegrationParameters(Points,SEC`MST,SEC`MST_Tau,N);
	vprint SE,2 : "#Abscissas:",#Abscissas;
	vprint SE,3 : "Abscissas:",Abscissas;
	vprint SE,2 : "#Weights:",#Weights;
	vprint SE,3 : "Weights:",Weights;
	vprint SE,2 : "StepLength:",StepLength;
	vprint SEP,1 : "Precision Abscissas:",Precision(Abscissas[1]);
	vprint SEP,1 : "Precision Weights:",Precision(Weights[1][1]);
	vprint SEP,1 : "Precision StepLength:",Precision(StepLength);
	Cputime(t);

	// Integrals
	vprint SE,1 : "Integrating...";
	t := Cputime();
	Integrals := [];
	for k in [1..d-1] do
		vprint SE,2 : "Integrating edge nr.",k;
		I := RS_SEIntegrate( SEC`MST[k],Points,DFF,Abscissas,Weights,StepLength,d,N );
		Append(~Integrals,I);
	end for;
	Cputime(t);

	// Matrix of periods
	//PM_ := ZeroMatrix(CC,g,2*g);
	PM_ := ZeroMatrix(CC,g,(d-1)*(N-1));
	//PM_ := ZeroMatrix(CC,g,2*g+d-1);
	for k in [1..d-1] do
		for l in [1..N-1] do
		//for l in [1..N] do
			Ind2 := (N-1)*(k-1) + l;
			//Ind2 := N*(k-1) + l;
			for j in [1..g] do
				//PM_[j][Ind2] := Zeta^(l-1) * (1 - Zeta^(N-1)) * Integrals[k][j];
				//PM_[j][Ind2] := Zeta^(l-1) * (1 - Zeta) * Integrals[k][j];
				PM_[j][Ind2] := Integrals[k][j][l];
			end for;
		end for;
	end for;
	vprint SE,3: "Integrals:",PM_;

	// Intersection matrix
	vprint SE,1 : "Computing skd-Matrix...";
	t := Cputime();
	skd_Matrix := [ [] : j in [1..d-1] ];
	for j in [1..d-1] do
		skd_Matrix[j][j] := <0,0,0>;
		for l in [j+1..d-1] do
			skd := RS_SEIntersection(SEC`MST[j],SEC`MST[l],Points,N);
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
	CF, SEC`SymplecticTransformation := RS_SymplecticBasis(SEC`IntersectionMatrix);
	Cputime(t);
	vprint SE,2: "ST:",SEC`SymplecticTransformation;
	vprint SE,3: "CF:",CF;
	
	
	ST_CC := ChangeRing(Transpose(SEC`SymplecticTransformation),CC);
	vprint SE,1 : "Matrix multiplication 1...";
	t := Cputime();
	PMAPMB := PM_ * ST_CC;
	Cputime(t);
	
	PM_A := ColumnSubmatrixRange(PMAPMB,1,g);
	PM_B := ColumnSubmatrixRange(PMAPMB,g+1,2*g);

	vprint SE,3 : "Dependent columns:",ColumnSubmatrixRange(PMAPMB,2*g+1,Nrows(ST_CC));
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
				if Abs(Re(PM[j][k])) lt RS_GetGlobalError() then
					PM[j][k] := Zero(C_) + Im(PM[j][k])*i;
				end if;
				if Abs(Im(PM[j][k])) lt RS_GetGlobalError() then
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
		//if MaxSymDiff ge RS_GetGlobalError() then
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
end intrinsic;
