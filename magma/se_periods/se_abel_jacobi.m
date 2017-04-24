/*******************************************************************************

 Copyright (C) 2017
 Christian Neurohr
 
*******************************************************************************/

// Import global settings
import "se_period_matrix.m": SE_PeriodMatrix;
import "se_de_int.m": DE_Integrals_Edge_AJM, DE_Int_Params_AJM;
import "se_help_funcs.m": MakeCCVector, Distance;
import "se_anal_cont.m": AC_mthRoot;

C20<i> := ComplexField(20);

function PeriodLatticeReduction(V,SEC)
// Reduce V \in \C^g modulo period matrix (A,B)
	if Universe(V) eq SEC`ComplexField then
		V := [ Re(v) : v in V ] cat [ Im(v) : v in V ];
	end if;
	// Assume now V in \R^2g
	W := ChangeRing(SEC`ReductionMatrix,SEC`RealField) * Matrix(SEC`RealField,2*SEC`Genus,1,V);
	return [w - Round(w) : w in Eltseq(W)];
end function;


procedure SE_ReductionMatrix(SEC)
// Computes a matrix for reduction modulo the period lattice given by the big matrix (A,B)
	if not assigned SEC`ReductionMatrix then
		PM_AB := SEC`BigPeriodMatrix;
		g := SEC`Genus;
		M := Matrix(SEC`RealField,2*g,2*g,[]);
		for j in [1..g] do
			for k in [1..g] do
				M[j,k] := Re(PM_AB[j,k]);
				M[j+g,k] := Im(PM_AB[j,k]);
				M[j,k+g] := Re(PM_AB[j,k+g]);	
				M[j+g,k+g] := Im(PM_AB[j,k+g]);
			end for;
		end for;
		// Matrix inversion
		SEC`ReductionMatrix := M^(-1);
	end if;
end procedure;

procedure SE_TreeMatrix(SEC)
// Computes a matrix describing paths in the spanning tree from P_0 => P_i for all ramification points P_i }
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


intrinsic SE_PrincipalBranch( x::FldComElt, SEC::SECurve : Fiber := false ) -> SeqEnum[FldComElt]
{ Returns y(x) = f(x)^(1/m) }
	m := SEC`Degree[1];
	fx := SEC`LeadingCoeff * &*[ (x-P) : P in SEC`BranchPoints ];
	if m eq 2 then
		y := Sqrt(fx);
	else
		y := Exp((1/m)*Log(fx));
	end if;
	if Fiber then
		return [ SEC`Zetas[2*k-1]*y : k in [1..m] ];
	else
		return y;
	end if;
end intrinsic;
	
				
procedure SE_RamificationPoints_AJM( SEC )
// Computes the Abel-Jacobi map of the divisors D = [ P_i - P_0 ] for i = 1,...,n
// and of D = [ \sum_{j = 0}^{delta} P_{\infty}^j - \delta P_0 ] 
	C<i> := SEC`ComplexField;
	Q := Rationals();
	g := SEC`Genus;
	m := SEC`Degree[1];
	n := SEC`Degree[2];
	if not assigned SEC`AJM_RamificationPoints then
		AJMRamPoints := [];
		for j in [1..n] do
			V := Matrix(C,g,1,[]);
			TreePath := SEC`TreeMatrix[j];
			for k in [1..n-1] do
				V +:= TreePath[k] * SEC`ElementaryIntegrals[k];
			end for;
			V_red := PeriodLatticeReduction(ElementToSequence(V),SEC);
			Append(~AJMRamPoints,Matrix(Q,2*g,1,[ BestApproximation(v,m) : v in V_red ]));
			
		end for;
		SEC`AJM_RamificationPoints := AJMRamPoints;
	end if;
	if not assigned SEC`AJM_SumOfInftyPoints then
		delta,a,b := Xgcd(m,n);
		SumOfInftyPoints := b * &+[ V : V in SEC`AJM_RamificationPoints ];
		SEC`AJM_SumOfInftyPoints := Matrix(Q,2*SEC`Genus,1,[ z - Round(z) : z in Eltseq(SumOfInftyPoints)] );
	end if;
end procedure;


intrinsic SE_AbelJacobi( D::SEDivisor, SEC::SECurve ) -> Mtrx
{ Computes the Abel-Jacobi map with basepoint P_0 of the zero divisor [D - deg(D)*P_0] }

	require D`Curve eq <SEC`Degree[1],SEC`DefiningPolynomial>: "Divisor has to be defined on points of the curve.";

	// Fields
	C<i> :=  SEC`ComplexField;
	R := SEC`RealField; 
	Q := RationalField();          
	
	// Split up computation
	ComplexIntegral := Matrix(C,SEC`Genus,1,[]);
	RationalIntegral := Matrix(Q,2*SEC`Genus,1,[]);
	RealIntegral := Matrix(R,2*SEC`Genus,1,[]);
			
	ComplexEdges := [];
	// Sort out infinite points and split up the rational part of the integrals
	for P in D`Support do
		v_P :=  P[2];
		if #P[1] eq 1 then
			// Infinite points
			if SEC`Degree[3] eq 1 then
				RationalIntegral +:= v_P * SEC`AJM_InftyPoints[1];
			else
				SE_InfinitePoints_AJM(Floor(P[1][1]),SEC);
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

	// Compute integration parameters
	DEInt := DE_Int_Params_AJM( ComplexEdges, SEC );
	
	// Actual integrations from P_k to P
	for CE in ComplexEdges do
		ComplexIntegral +:= CE[2] * Matrix(C,SEC`Genus,1,DE_Integrals_Edge_AJM(CE,SEC,DEInt));
	end for;

	// Differential change matrix
	ComplexIntegral := SEC`DifferentialChangeMatrix * ComplexIntegral;

	// Reduce complex vector modulo (A,B)
	TotalRealIntegral := Matrix(R,2*SEC`Genus,1,PeriodLatticeReduction(Eltseq(ComplexIntegral),SEC));

	// Add integrals corresponding to ramification points and infinite points 
	TotalRealIntegral +:= ChangeRing(RationalIntegral,R) + RealIntegral;
			
	// Reduce again modulo \Z and return results
	//R_ := RealField(SEC`Prec);
	return Matrix(R,2*SEC`Genus,1,[v - Round(v) : v in Eltseq(TotalRealIntegral) ]);
	
end intrinsic;





		
// Integrating to infty // not working properly // use moving lemma instead

function AC_mthRootInfty(x,p,Points,Zeta,N:long:=false)
// Analytic continuation of y = f(x)^(1/N)  for x -> \infty on one sheet, avoiding branch cuts
	long := true;
	if #Points eq 1 or long then
		prod := 1;
		for k in [1..#Points] do
			prod *:= ( 2*p + (x-1) * Points[k]  )^(1/N);
		end for;
		return prod;
	end if;
end function;

/*intrinsic SE_AJInfinity( SEC::SECurve : Ind := 1 ) -> SeqEnum[FldComElt]
{ Integrate from P_0 to \infty }

	"################## Integrate P_0 to Infinity ####################";

	C<i> := SEC`ComplexField; C_0 := Zero(C);
	Points := SEC`BranchPoints;
	m := SEC`Degree[1];
	n := SEC`Degree[2];
	g := SEC`Genus;
	print "Points:",Points;

	// Basepoint
	x_0 := Points[SEC`Basepoint];
	print "x_0:",x_0;

 	MaxIm := Max( [Abs(Im(Points[k])) : k in [1..n] ]);
	print "MaxIm:",MaxIm;

	x_P := 2*(Min(Re(x_0)-1,-1) + i * (MaxIm + 1));
	//p := CC!(1.6666666666666666667 - 1.6665612762750131730*i);
	//p := Re(x_0) + 1000*i*(MaxIm+1);

	InfTau := RS_SEInfTau( x_P, Points : Lambda := Pi(C20)/2 );
	print "InfTau:",InfTau;
	// TODO Better Tau?
	if InfTau lt SEC`Tau then
		// Update integration parameters
		print "Old tau:",SEC`Tau;
		SEC`Tau := InfTau;
		print "New tau:",SEC`Tau;
		SEC`Abscissas, SEC`Weights, SEC`StepLength := RS_SEDEIntegrationParameters(Points,SEC`MST,InfTau,m,SEC`ComplexField);
	end if;

	Y_P := RS_PrincipalBranch(SEC,x_P:Fiber);
	
	P := [x_P,Y_P[Ind]];
	print "P:",P;

	// Integrate from P_0 to P
	V1 := RS_AJMIntegrate(SEC,1,P);

	// Integrate from P to \infty
	Fact := 2*x_P;
	print "Fact:",Fact;

	// Vector integral
	Integral := [ C_0 : j in [1..g] ];

	// Differentials
	DFF := &cat[ DFF_i : DFF_i in SEC`HolomorphicDifferentials ];

	// Initiate on x = 0, dx = 1
	z := AC_mthRootInfty(0,x_P,Points,SEC`Zetas,m);
	for j in [1..g] do
		w := DFF[j];
		Denom := z^w[2];
		Integral[j] +:= (1/Denom);
	end for;

	for t in [1..#SEC`Abscissas] do
		x := SEC`Abscissas[t];

		// Analytic continuation
		z1 := AC_mthRootInfty(x,x_P,Points,SEC`Zetas,m);
		z2 := AC_mthRootInfty(-x,x_P,Points,SEC`Zetas,m);
		
		Enum1 := (1 - x);
		Enum2 := (1 + x);

		for j in [1..g] do
			w := DFF[j];
			Pow := n*w[2]/m - w[1] - 1;
			assert Pow ge 0;
			Denom1 := z1^w[2];
			Denom2 := z2^w[2];
			dx := SEC`Weights[t][1];
			//dx := SEC`Weights[t][1] * SEC`Weights[t][2]^2;
			//Integral[j] +:= ( Enum1^w[1] * c1^w[2] / Denom1 + Enum2^w[1] * c2^w[2] / Denom2) * dx;
			Integral[j] +:= ( (Enum1^Pow * Enum2 / Denom1) + (Enum2^Pow * Enum1 / Denom2) ) * dx;
		end for;
	end for;


	P_y_AC := AC_mthRootInfty(-1,x_P,Points,SEC`Zetas,m); // * 2^(d/N)
	//print "P_y_AC:",P_y_AC;
	print "Arg(P_y_AC):",Arg(P_y_AC);
	//print "P_y:",P[2];
	print "Arg(P_y):",Arg(P[2]);	

	k_ := (m/(2*SEC`Pi)) * ( Arg(P[2]) - Arg(P_y_AC) - Arg(SEC`LeadingCoeff) );
	print "k_:",k_;
	k := Round(k_);
	assert Abs(k - k_) lt 10^-10; // Check: k \in \Z ?
	for j in [1..g] do
		w := DFF[j];
		Pow := -2*w[2]*k mod (2*m) + 1;
		print "Pow:",Pow;
		Integral[j] *:= SEC`LeadingCoeff^w[2] * SEC`Zetas[Pow] * SEC`StepLength * Fact^(w[1]+1);
	end for;

	print "Integral_P0toP",V1;
	print "Integral_PtoInfty:",Integral;

	V := ElementToSequence(V1);
	// I_{\infty} = (I_{P_0,P} + I_{P,\infty})
	Integral_P0toInf := [ V[j] + Integral[j] : j in [1..g] ];

	print "Integral_P0toInfty:",Integral_P0toInf;

	M := Matrix(C,g,1,Integral_P0toInf);
	print "m*M:",m*M;
	RedM := RS_LatticeReduction(SEC,ElementToSequence(m*M));
	print "RedM:",RedM;

	ReducedIntegral := RS_LatticeReduction(SEC,Integral_P0toInf);

	print "ReducedIntegral:",ReducedIntegral;

	//Res := [ BestApproximation(ReducedIntegral[j][1],m) : j in [1..2*g] ];
	Res := [ ReducedIntegral[j][1] : j in [1..2*g] ];
	return Res;
end intrinsic;*/
