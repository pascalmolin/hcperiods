/*******************************************************************************

 Copyright (C) 2017 
 Code adapted from Pascal Molin,
 written by Christian Neurohr

*******************************************************************************/

// Import functions;
import "se_spanning_tree.m": SpanningTree;
import "se_de_int.m": DE_Int_Params, DE_Integrals_Tree;
import "se_gc_int.m": GC_Integrals_Tree;
import "se_anal_cont.m:": AC_mthRoot;
import "se_intersection.m": SE_IntersectionMatrix, SE_IntersectionNumber;
import "se_symplectic_basis.m":SymplecticBasis;


intrinsic SE_BigPeriodMatrix( f::RngMPolElt : Prec := 40 ) -> Mtrx
{ Computes a big period matrix associated to f(x,y) = y^m - g(x) to precision Prec }
	S := SE_Curve(f:Prec:=Prec,Small:=false,AbelJacobi:=false);
	return S`BigPeriodMatrix;
end intrinsic;
intrinsic SE_BigPeriodMatrix( f::RngUPolElt,m::RngIntElt : Prec := 40 ) -> Mtrx
{ Computes a big period matrix associated to y^m = f(x) to precision Prec }
	S := SE_Curve(f,m:Prec:=Prec,Small:=false,AbelJacobi:=false);
	return S`BigPeriodMatrix;
end intrinsic;
intrinsic SE_BigPeriodMatrix( Points::SeqEnum[FldComElt],m::RngIntElt : LeadingCoeff := 1, Prec := 40 ) -> Mtrx
{ Computes a big period matrix associated to y^m = LeadingCoeff * \prod[(x-p) : p in Points] to precision Prec }
	S := SE_Curve(Points,m:Prec:=Prec,Small:=false,AbelJacobi:=false);
	return S`BigPeriodMatrix;
end intrinsic;
intrinsic SE_SmallPeriodMatrix( f::RngMPolElt : Prec := 40 ) -> Mtrx
{ Computes a small period matrix associated to f(x,y) = y^m - g(x) to precision Prec }
	S := SE_Curve(f:Prec:=Prec,Small:=true,AbelJacobi:=false);
	return S`SmallPeriodMatrix;
end intrinsic;
intrinsic SE_SmallPeriodMatrix( f::RngUPolElt,m::RngIntElt : Prec := 40 ) -> Mtrx
{ Computes a small period matrix associated to y^m = f(x) to precision Prec }
	S := SE_Curve(f,m:Prec:=Prec,Small:=true,AbelJacobi:=false);
	return S`SmallPeriodMatrix;
end intrinsic;
intrinsic SE_SmallPeriodMatrix( Points::SeqEnum[FldComElt],m::RngIntElt : LeadingCoeff := 1, Prec := 40 ) -> Mtrx
{ Computes a small period matrix associated to y^m = LeadingCoeff * \prod[(x-p) : p in Points] to precision Prec }
	S := SE_Curve(Points,m:Prec:=Prec,Small:=true,AbelJacobi:=false);
	return S`SmallPeriodMatrix;
end intrinsic;


procedure SE_PeriodMatrix( SEC : Small := true, ReductionMatrix := false )
// Computes period matrices associated to the SE_Curve object SEC }

	if &and[Small,not assigned SEC`SmallPeriodMatrix] or &and[not Small, not assigned SEC`BigPeriodMatrix] then 

	// Degrees
	m := SEC`Degree[1]; n := SEC`Degree[2];

	// Complex fields
	C<i> := SEC`ComplexField;
	CS<i> := ComplexField(SEC`Prec);
	
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
		DEInts := DE_Integration(STree`Params,SEC);
		vprint SE,1 : "using double-exponential integration...";
		vprint SE,2 : [ <D`NPoints,D`r> : D in DEInts ];	
		Periods, ElemInts := DE_Integrals_Tree(SEC,DEInts);
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
	BPM := SEC`DifferentialChangeMatrix * PM_AB;
	// Compute reduction matrix
	if ReductionMatrix then
		// Embed big period matrix in \R^{2g \times 2g}
		M := Matrix(Parent(Re(BPM[1][1])),2*g,2*g,[]);
		for j in [1..g] do
			for k in [1..g] do
				M[j,k] := Re(BPM[j,k]);
				M[j+g,k] := Im(BPM[j,k]);
				M[j,k+g] := Re(BPM[j,k+g]);	
				M[j+g,k+g] := Im(BPM[j,k+g]);
			end for;
		end for;
		// Matrix inversion
		SEC`ReductionMatrix := ChangeRing(M^(-1),SEC`RealField);
	end if;

	if Small then
		vprint SE,1 : "Matrix inversion...";
		PM_AInv := ColumnSubmatrixRange(BPM,1,g)^(-1);
		vprint SE,1 : "Matrix multiplication 2...";
		Tau := PM_AInv * ColumnSubmatrixRange(BPM,g+1,2*g);
		Tau := ChangeRing(Tau,CS);

		// Assign small period matrix
		SEC`SmallPeriodMatrix :=  Tau;

		// Testing for symmetry of the period matrix
		vprint SE,1 : "Testing symmetry...";
		MaxSymDiff := 0;
		for j in [1..g] do
			for k in [j+1..g] do
				MaxSymDiff := Max(MaxSymDiff,Abs(Tau[j][k] - Tau[k][j]));
			end for;
		end for;
		vprint SE,2 : "Maximal symmetry deviation:",MaxSymDiff;
		if MaxSymDiff ge SEC`Error then
			print "Small period matrix: Requested accuracy could not not be reached.";
			print "Significant digits:",Floor(-Log(10,MaxSymDiff));
		end if;	
		
		// Testing positive definiteness of the imaginary part of the period matrix
		if false then
			vprint SE,1 : "Testing positive definite...";
			Tau_Im := ZeroMatrix(SEC`RealField,g,g);
			for j in [1..g] do
				for k in [j..g] do
					Tau_Im[j][k] := Real(Im(Tau[j][k]));
					Tau_Im[k][j] := Tau_Im[j][k];
				end for;
			end for;
			assert IsPositiveDefinite(Tau_Im);
		end if;
		
	end if;
	
	// Assign big period matrix
	SEC`BigPeriodMatrix := ChangeRing(BPM,CS);

	end if;
end procedure;
