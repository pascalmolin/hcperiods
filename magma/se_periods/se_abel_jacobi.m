/*******************************************************************************

 Copyright (C) 2017
 Christian Neurohr
 
*******************************************************************************/

// Import global settings
import "se_period_matrix.m": SE_PeriodMatrix;
import "se_de_int.m": DE_Integrals_Edge_AJM, DE_Int_Params;
import "se_help_funcs.m": MakeCCVector, Distance;
import "se_anal_cont.m": AC_mthRoot;
import "se_spanning_tree.m": DE_Params_AJM;


function PeriodLatticeReduction(V,SEC)
// Reduce V \in \C^g modulo period matrix (A,B)
	if Universe(V) eq SEC`ComplexField then
		V := [ Re(v) : v in V ] cat [ Im(v) : v in V ];
	end if;
	// Assume now V in \R^2g
	W := SEC`ReductionMatrix * Matrix(BaseRing(Parent(SEC`ReductionMatrix)),2*SEC`Genus,1,V);
	return [w - Round(w) : w in Eltseq(W)];
end function;


procedure SE_ReductionMatrix(SEC)
// Computes a matrix for reduction modulo the period lattice given by the big period matrix (A,B)
	if not assigned SEC`ReductionMatrix then
		BPM := SEC`BigPeriodMatrix;
		g := SEC`Genus;
		// Embed big period matrix in R^{2g \times 2g}
		M := Matrix(SEC`RealField,2*g,2*g,[]);
		for j in [1..g] do
			for k in [1..g] do
				M[j,k] := Re(BPM[j,k]);
				M[j+g,k] := Im(BPM[j,k]);
				M[j,k+g] := Re(BPM[j,k+g]);	
				M[j+g,k+g] := Im(BPM[j,k+g]);
			end for;
		end for;
		// Matrix inversion
		SEC`ReductionMatrix := M^(-1);
	end if;
end procedure;


procedure SE_TreeMatrix(SEC)
// Computes a 'tree-matrix' which contains paths in the spanning tree from P_0 => P_i for all ramification points P_i
	if not assigned SEC`TreeMatrix then
	n := SEC`Degree[2];
	TM := ZeroMatrix(Integers(),n,n-1);
	Taken := [ 0 : j in [1..n] ];
	Tree := SEC`SpanningTree`Edges;
	P_0 := Tree[1][1];
	Taken[P_0] := 1;
	for j in [1..n-1] do
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
	// Shift by basepoint for Abel-Jacobi map
	P_0 := SEC`Basepoint; 
	P_0P_0 := TM[P_0];
	for j in [1..n] do
		TM[j] -:= P_0P_0;
	end for;
	SEC`TreeMatrix := TM;
	end if;
end procedure;


intrinsic SE_Branches( x0::FldComElt, SEC::SECurve : Principal := false ) -> SeqEnum[FldComElt]
{ Computes the solutions of y^m = f(x0) }
	if Distance(x0,SEC`BranchPoints) lt SEC`Error then
		y := 0;
	else
		m := SEC`Degree[1];
		fx0 := SEC`LeadingCoeff * &*[ (x0-P) : P in SEC`BranchPoints ];
		if m eq 2 then
			y := Sqrt(fx0);
		else
			y := Exp((1/m)*Log(fx0));
		end if;
	end if;
	if Principal then
		return y;
	else
		return [ SEC`Zetas[2*k-1]*y : k in [1..m] ];
	end if;
end intrinsic;
	
				
procedure SE_RamificationPoints_AJM( SEC )
// Computes the Abel-Jacobi map of the divisors D_i = [ P_i - P_0 ] for i = 1,...,n
// and of D = [ \sum_{j = 0}^{delta} P_{\infty}^j - \delta P_0 ] 
	m := SEC`Degree[1];
	n := SEC`Degree[2];
	if not assigned SEC`AJM_RamificationPoints then
		AJMRamPoints := [];
		for j in [1..n] do
			V := Matrix(SEC`ComplexField,SEC`Genus,1,[]);
			TreePath := SEC`TreeMatrix[j];
			for k in [1..n-1] do
				V +:= TreePath[k] * SEC`ElementaryIntegrals[k];
			end for;
			V_red := PeriodLatticeReduction(ElementToSequence(V),SEC);
			Append(~AJMRamPoints,Matrix(Rationals(),2*SEC`Genus,1,[ BestApproximation(v,m) : v in V_red ]));
			
		end for;
		SEC`AJM_RamificationPoints := AJMRamPoints;
	end if;
	if not assigned SEC`AJM_SumOfInftyPoints then
		delta,a,b := Xgcd(m,n);
		SumOfInftyPoints := b * &+[ V : V in SEC`AJM_RamificationPoints ];
		SEC`AJM_SumOfInftyPoints := Matrix(Rationals(),2*SEC`Genus,1,[ z - Round(z) : z in Eltseq(SumOfInftyPoints)] );
	end if;
end procedure; 


intrinsic SE_AbelJacobi( P::SeqEnum, Q::SeqEnum, SEC::SECurve ) -> Mtrx
{ Computes the Abel-Jacobi map of the zero divisor [ Q - P ] }
	Points := [<Q,1>, <P,-1>];
	return SE_AbelJacobi(SE_Divisor(Points,SEC),SEC);
end intrinsic;


intrinsic SE_AbelJacobi( D::SEDivisor, P0::SeqEnum, SEC::SECurve ) -> Mtrx
{ Computes the Abel-Jacobi map of the divisor D with basepoint P0 }
	if D`Degree ne 0 then 
		D0 := SE_Divisor([<P0,-D`Degree>],SEC);
		return SE_AbelJacobi(D+D0,SEC);
	else
		return SE_AbelJacobi(D,SEC);
	end if;
	
end intrinsic;


intrinsic SE_AbelJacobi( D::SEDivisor, SEC::SECurve ) -> Mtrx
{ Computes the Abel-Jacobi map of the zero divisor D }

	require D`Curve eq <SEC`Degree[1],SEC`DefiningPolynomial>: "Divisor has to be defined on points of the curve.";
	require D`Degree eq 0 : "Degree of divisor has to be zero.";

	// Fields
	C<I> :=  SEC`ComplexField;
	R := SEC`RealField; 
	Q := RationalField();          
	
	// Split up computation
	RationalIntegral := Matrix(Q,2*SEC`Genus,1,[]);
	RealIntegral := Matrix(R,2*SEC`Genus,1,[]);
			
	ComplexEdges := [];
	// Sort out infinite points and split up the rational part of the integrals
	for P in D`Support do
		v_P :=  P[2];
		if #P[1] eq 1 then
			// Infinite points
			SE_AJM_InftyPoints(Floor(P[1][1]),SEC);
			if SEC`Degree[3] eq 1 then
				RationalIntegral +:= v_P * SEC`AJM_InftyPoints[1];
			else
				RealIntegral +:= v_P * SEC`AJM_InftyPoints[Floor(P[1][1])];
			end if;
		else
			// Finite points
			Dist, Ind := Distance(P[1][1],SEC`BranchPoints);  
			RationalIntegral +:= v_P * SEC`AJM_RamificationPoints[Ind];
			if Dist gt SEC`Error then
				Append(~ComplexEdges,Append(P,Ind));
			end if;
		end if;
	end for;

	// No integration needed in this special case
	if #ComplexEdges eq 0 then
		RealIntegral +:= ChangeRing(RationalIntegral,R);
		return Matrix(R,2*SEC`Genus,1,[v - Round(v) : v in Eltseq(RealIntegral) ]);
	end if;

	// Integration parameters
	Params := DE_Params_AJM(ComplexEdges,SEC);
	vprint SE,2 : "Parameter(AJM):",Params;
	ExtraPrec := 2*Ceiling(Log(10,Params[1]/SEC`SpanningTree`Params[1]));
	if ExtraPrec gt 0 then
		vprint SE,2 :"Extra precision (AJM):",ExtraPrec;
		C<I> := ComplexField(Precision(C)+ExtraPrec);
		SEC`ComplexField := C;
		// Update curve data
		f_x := ChangeRing(SEC`DefiningPolynomial,C);
		SEC`ComplexPolynomial := f_x;
		// Roots of f_x
		Roots_fx := Roots(f_x);
		SEC`BranchPoints := [ R[1] : R in Roots_fx ];
		//SEC`DifferentialChangeMatrix := ChangeRing(SEC`DifferentialChangeMatrix,C);
	end if;
	DEInt := DE_Int_Params(Params,SEC:AJM:=true);

	// Actual integrations from P_k to P
	ComplexIntegral := Matrix(C,SEC`Genus,1,[]);
	for CE in ComplexEdges do
		ComplexIntegral +:= CE[2] * Matrix(C,SEC`Genus,1,DE_Integrals_Edge_AJM(CE,SEC,DEInt));
	end for;

	// Differential change matrix
	ComplexIntegral := ChangeRing(SEC`DifferentialChangeMatrix,C) * ComplexIntegral;

	// Reduce complex vector modulo (A,B)
	TotalRealIntegral := Matrix(R,2*SEC`Genus,1,PeriodLatticeReduction(Eltseq(ComplexIntegral),SEC));

	// Add integrals corresponding to ramification points and infinite points 
	TotalRealIntegral +:= ChangeRing(RationalIntegral,R) + RealIntegral;
			
	// Reduce again modulo \Z and return results
	return Matrix(R,2*SEC`Genus,1,[v - Round(v) : v in Eltseq(TotalRealIntegral) ]);
	
end intrinsic;
