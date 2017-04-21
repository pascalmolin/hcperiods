/*******************************************************************************

 Copyright (C) 2017 
 Code adapted from Pascal Molin,
 written by Christian Neurohr

*******************************************************************************/

// Import functions;
import "se_spanning_tree.m": SpanningTree;
import "se_de_int.m": DE_Int_Params_Tree, DE_Integrals_Tree;
import "se_gc_int.m": GC_Integrals_Tree;
import "se_anal_cont.m:": AC_mthRoot;
import "se_intersection.m": SE_IntersectionMatrix, SE_IntersectionNumber;
import "se_symplectic_basis.m":SymplecticBasis;


// Verbose
declare verbose SE,3;

intrinsic SE_BigPeriodMatrix( f::RngMPolElt : Prec := 20 ) -> Mtrx
{ Computes a big periodmatrix associated to f(x,y) = y^m - g(x) to precision Prec }
	S := SE_Curve(f:Prec:=Prec,Small:=false,AbelJacobi:=false);
	return S`BigPeriodMatrix;
end intrinsic;
intrinsic SE_BigPeriodMatrix( f::RngUPolElt,m::RngIntElt : Prec := 20 ) -> Mtrx
{ Computes a big periodmatrix associated to y^m = f(x) to precision Prec }
	S := SE_Curve(f,m:Prec:=Prec,Small:=false,AbelJacobi:=false);
	return S`BigPeriodMatrix;
end intrinsic;
intrinsic SE_SmallPeriodMatrix( f::RngMPolElt : Prec := 20 ) -> Mtrx
{ Computes a big periodmatrix associated to f(x,y) = y^m - g(x) to precision Prec }
	S := SE_Curve(f:Prec:=Prec,Small:=true,AbelJacobi:=false);
	return S`SmallPeriodMatrix;
end intrinsic;
intrinsic SE_SmallPeriodMatrix( f::RngUPolElt,m::RngIntElt : Prec := 20 ) -> Mtrx
{ Computes a big periodmatrix associated to y^m = f(x) to precision Prec }
	S := SE_Curve(f,m:Prec:=Prec,Small:=true,AbelJacobi:=false);
	return S`SmallPeriodMatrix;
end intrinsic;

procedure SE_PeriodMatrix( SEC : Small := true, AbelJacobi := true )
// Computes period matrices associated to the SE_Curve object SEC }

	if &and[Small,not assigned SEC`SmallPeriodMatrix] or &and[not Small, not assigned SEC`BigPeriodMatrix] then 

	// Degrees
	m := SEC`Degree[1]; n := SEC`Degree[2];

	// Complex fields and constants
	C20<i> := ComplexField(20);
	C<i> := SEC`ComplexField;
	CS<i> := ComplexField(SEC`Prec);
	C_0 := Zero(C);
	
	// Spanning tree
	STree := SEC`SpanningTree;
	vprint SE,1 : "Spanning tree:",STree;	

	// Holomorphic differentials
	DFF := SEC`HolomorphicDifferentials;
	vprint SE,1 : "Holomorphic differentials:",DFF;

	// Genus
	g := SEC`Genus;
	vprint SE,1 : "Genus:",g;

	// Computing periods
	vprint SE,1 : "Integrating...";
	if SEC`IntegrationType eq "DE" then
		// Integration parameters
		vprint SE,1 : "Computing integration parameters...";
		DEInt := DE_Int_Params_Tree(SEC);
		vprint SE,1 : "using double-exponential integration...";
		vprint SE,2 : DEInt;	
		Periods, ElemInts := DE_Integrals_Tree(SEC,DEInt);
	elif SEC`IntegrationType eq "GC" then
		vprint SE,1 : "using Gauss-Chebychev integration...";
		Periods, ElemInts := GC_Integrals_Tree(SEC);
	else
		error Error("Invalid integration type.");
	end if;

	// Elementary integrals
	SEC`ElementaryIntegrals := [];
	for k in [1..n-1] do
		V := Matrix(C,g,1,ElemInts[k]);
		Append(~SEC`ElementaryIntegrals,SEC`DifferentialChangeMatrix*V);
	end for;

	// Period matrix
	PM := ZeroMatrix(C,g,(m-1)*(n-1));
	for k in [1..n-1] do
		for l in [1..m-1] do
			Ind := (m-1)*(k-1) + l;
			for j in [1..g] do
				PM[j][Ind] := Periods[k][j][l];
			end for;
		end for;
	end for;

	// Intersection matrix
	vprint SE,1 : "Computing spsm-matrix...";
	spsm_Matrix := [ [] : j in [1..n-1] ];
	for j in [1..n-1] do
		spsm_Matrix[j][j] := <1,m-1>;
		for k in [j+1..n-1] do
			spsm := SE_IntersectionNumber(SEC`SpanningTree`Edges[j],SEC`SpanningTree`Edges[k],SEC`BranchPoints,m,n,SEC`Zetas);
			spsm_Matrix[j][k] := <spsm[1] mod m,spsm[2] mod m>;
		end for;
	end for;
	vprint SE,3: "spsm_Matrix:",spsm_Matrix;
	

	vprint SE,1 : "Computing intersection matrix...";
	SEC`IntersectionMatrix := SE_IntersectionMatrix(spsm_Matrix,m,n);
	assert Rank(SEC`IntersectionMatrix) eq 2*g;

	// Symplectic reduction of intersection matrix
	vprint SE,1 : "Performing symplectic reduction...";
	CF, ST := SymplecticBasis(SEC`IntersectionMatrix);
	ST_C := ChangeRing(Transpose(ST),C);
	vprint SE,3: "ST:",ST;
	vprint SE,3: "CF:",CF;

	vprint SE,1 : "Matrix multiplication 1...";
	PMAPMB := PM * ST_C;

	PM_AB := ColumnSubmatrixRange(PMAPMB,1,2*g);
	vprint SE,3 : "Dependent columns:",ColumnSubmatrixRange(PMAPMB,2*g+1,Nrows(ST));

	// Compute big period matrix
	BigPM := SEC`DifferentialChangeMatrix * PM_AB;
	SEC`BigPeriodMatrix := ChangeRing(BigPM,CS);

	if Small then
		vprint SE,1 : "Matrix inversion...";
		PM_AInv := ColumnSubmatrixRange(BigPM,1,g)^(-1);
		vprint SE,1 : "Matrix multiplication 2...";
		Tau := PM_AInv * ColumnSubmatrixRange(BigPM,g+1,2*g);
		Tau := ChangeRing(Tau,CS);

		// Testing for symmetry of the period matrix
		vprint SE,1 : "Testing symmetry...";
		MaxSymDiff := 0;
		for j in [1..g] do
			for k in [j+1..g] do
				MaxSymDiff := Max(MaxSymDiff,Abs(Tau[j][k] - Tau[k][j]));
			end for;
		end for;
		vprint SE,1 : "Maximal symmetry deviation:",MaxSymDiff;


		if MaxSymDiff ge SEC`Error then
			print "Warning: requested accuracy could not not be reached.";
			print "Significant digits:",Floor(-Log(10,MaxSymDiff));
		end if;	
		
		// Testing positive definiteness of the imaginary part of the period matrix
		vprint SE,1 : "Testing positive definite...";
		Tau_Im := ZeroMatrix(SEC`RealField,g,g);
		for j in [1..g] do
			for k in [j..g] do
				Tau_Im[j][k] := Real(Im(Tau[j][k]));
				Tau_Im[k][j] := Tau_Im[j][k];
			end for;
		end for;
		assert IsPositiveDefinite(Tau_Im);
		SEC`SmallPeriodMatrix := Tau;
	end if;

	end if;
end procedure;
