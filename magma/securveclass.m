// Import global settings
import "superelliptic.m": RS_SEIntegrate, RS_SEIntegrate3, RS_NRootAC, RS_SEMST, RS_SEMST2, RS_MakeCCVectors, RS_SETau3, RS_SEInfTau;
import "comparefunctions.m": RS_CompareByFirstComplexEntry,RS_CompareFldComElt;
import "globalprecision.m": RS_Config;
import "pathmethods.m": RS_ChoiceOfPathStar;
import "cvector.m": RS_Quotrem_I, RS_Quotrem_II;

C_20<i> := ComplexField(20);
C_Pi := Real(Pi(C_20));

// Superelliptic curves
declare type SECurve;

declare attributes SECurve: DefiningPolynomial, Genus, Degree, BranchPoints, LowPrecBranchPoints, MST,  Tau, ComplexField, RealField, HolomorphicDifferentials, Prec, Theta, SmallPeriodMatrix, BigPeriodMatrix, Pi, PeriodMatrix, ReductionMatrix, Abscissas, Weights, StepLength, AJM_RamPoints, AJM_InftyPoints, AbelJacobi, TreeMatrix, Error,
ElementaryIntegrals, spsm_Matrix, IntersectionMatrix, Basepoint, Zetas, LeadingCoeff;

// Constructor
intrinsic RS_SECurve( f::RngMPolElt : Prec := -1 ) -> SECurve
{ Creates the Riemann surface defined by f(x,y) = 0 and computes its function field and genus }
	
	// Requirements
	m := Degree(f,2); n := Degree(f,1); delta := Gcd(m,n);
	require &and[m ge 2, n ge 3] : "Degrees not supported."; 

	// Create object
	SEC := New(SECurve);	
	
	// Defining polynomial
	SEC`DefiningPolynomial := f;

	// Precision
	if Prec lt 20 then
		SEC`Prec := 20;
	else
		SEC`Prec := Prec;
	end if;

	// Error
	SEC`Error := 10^-(SEC`Prec+1);

	// Degrees
	SEC`Degree := [m,n,delta];	

	// Genus
	SEC`Genus := Round((1/2) * ((m-1)*(n-1) - delta + 1));

	// Complex field
	CC<i> := ComplexField(SEC`Prec+10);
	SEC`ComplexField := CC; SEC`Pi := Real(Pi(CC));

	// Real field
	SEC`RealField := RealField(SEC`Prec+5);
	
	// Univariate Polynomial
	f_x := ChangeRing(UnivariatePolynomial(Evaluate(f,[Parent(f).1,0])),CC);

	// Leading coefficient
	SEC`LeadingCoeff := Exp( (1/m) * -Log(-LeadingCoefficient(f_x)));

	// Branch Points
	SEC`BranchPoints := RS_Roots(f_x);

	// Low precision branch Points
	SEC`LowPrecBranchPoints := [ C_20!x : x in SEC`BranchPoints ];

	// Holomorphic differentials and genus
	SEC`HolomorphicDifferentials := RS_SEDifferentials(m,n:Int3:=false);

	// Root of unity powers
	SEC`Zetas := [ Exp(k*SEC`Pi*i/m) : k in [0..2*m-1] ];

	// Choose branch point with smallest real part as base point
	SEC`Basepoint := 1;

	// Compute everything
	RS_InitiateAbelJacobi(SEC);

	return SEC;
end intrinsic;

// Printing
intrinsic Print( SEC::SECurve )
{ Print Riemann surface }
	print "Superelliptic curve of genus", SEC`Genus ,"defined by: 0 =",SEC`DefiningPolynomial,"and prescribed precision",SEC`Prec;
	print "";
	print "Computed data:";
	print " Complex field:", SEC`ComplexField;
	print " Real field:", SEC`RealField;
	print " Branch points:", assigned SEC`BranchPoints;
	print " Roots of unity:", assigned SEC`Zetas;
	print " Maximal spanning tree:", assigned SEC`MST;
	print " Tau for DE-integration:", assigned SEC`Tau;
	print " spsm_Matrix:", assigned SEC`spsm_Matrix;
	print " Intersection matrix:", assigned SEC`IntersectionMatrix;
	print " Reduction matrix:", assigned SEC`ReductionMatrix;
	print " Basis of holomorphic differentials:", assigned SEC`HolomorphicDifferentials;
	print " Period matrix C:", assigned SEC`PeriodMatrix;
	print " Big period matrix (A,B):",assigned SEC`BigPeriodMatrix;
	print " Small period matrix A^(-1)B:",assigned SEC`SmallPeriodMatrix;
	print " Theta function:",assigned SEC`Theta;
	print " Abel-Jacobi map:", assigned SEC`AbelJacobi;
	print " AJM between ramification points:", assigned SEC`AJM_RamPoints;
	print " AJM to infinity:", assigned SEC`AJM_InftyPoints;
	print " Basepoint:", assigned SEC`Basepoint;
	print " Elementary integrals:", assigned SEC`ElementaryIntegrals;
	print " Tree matrix:", assigned SEC`TreeMatrix;
	print " Reduction matrix:", assigned SEC`ReductionMatrix;
end intrinsic;


function RS_AJMTau(x,j,Points)
	Tau := 4.0;
	for k in [1..#Points] do
		if k ne j then
			Tau := Min(Tau,RS_SETau3(x,Points[j],Points[k]));
		end if;
	end for;
	return Tau;
end function;
/*
function RS_FindPaths( SEC, D, NormalPoints : ED := true )
// Choose path such that tau for integration is maximal or euclidean distance is minimal
	n := SEC`Degree[2]; LowPrecBranchPoints := [ C_20!x : x in SEC`BranchPoints ]; x_P := C_20!x_P; BestTau := 0; ClosestBranchPoints := [];
	if ED then
		for j in NormalPoints do
			Dist, BestInd := RS_Distance(C_20!P[1],LowPrecBranchPoints);
			print "Dist:",Dist;
			print "Ind:",BestInd;

		BestTau := 4.0; // > Pi/2
		for k in [1..n] do
			if BestInd ne k then
				BestTau := Min(BestTau,RS_SETau3(x_P,Points[BestInd],Points[k]));
			end if;
		end for;
	else
		for j in [1..n] do
			Tau_k := 4.0; // > Pi/2
			for k in [1..n] do
				if k ne j then
					Tau_k := Min(Tau_k,RS_SETau3(x_P,Points[j],Points[k]));
				end if;
			end for;
			if Tau_k gt BestTau then
				BestInd := j;
				BestTau := Tau_k;
			end if;
		end for;
	end if;
	if BestTau lt SEC`Tau then
		// Update integration parameters
		SEC`Tau := BestTau;
		SEC`Abscissas, SEC`Weights, SEC`StepLength := RS_SEDEIntegrationParameters(Points,SEC`MST,BestTau,SEC`Degree[1],SEC`ComplexField);
	end if;
	return BestInd;	
end function;*/

intrinsic RS_LatticeReduction( SEC::SECurve, V::SeqEnum[FldComElt] ) -> Mtrx
{ Wrapper }
	assert #V eq SEC`Genus;
	return RS_LatticeReduction(SEC,[ Re(v) : v in V ] cat [ Im(v) : v in V ]);
end intrinsic;
intrinsic RS_LatticeReduction( SEC::SECurve, V::SeqEnum[FldReElt] ) -> SeqEnum[FldReElt]
{ Reduce V \in \C^g modulo period matrix (A B) }
	R := SEC`RealField; g := SEC`Genus;
	// Assume V in \R^2g
	W := SEC`ReductionMatrix * Matrix(R,2*g,1,V);
	return [w - Round(w) : w in ElementToSequence(W)];
end intrinsic;

procedure RS_ReductionMatrix( SEC : Recompute := false )
// Computes a matrix for reduction modulo big period matrix (A,B)
	if not assigned SEC`ReductionMatrix or Recompute then
		if not assigned SEC`BigPeriodMatrix then
			RS_SEPM(SEC:Small:=false,Recompute:=Recompute);
		end if;
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

intrinsic RS_PrincipalBranch( SEC::SECurve, x::FldComElt : Fiber := false ) -> SeqEnum[FldComElt]
{ Returns y(x) = f(x)^(1/m) }
	y := (1/SEC`LeadingCoeff) * (&*[ x-SEC`BranchPoints[k] : k in [1..SEC`Degree[2]] ])^(1/SEC`Degree[1]);
	if Fiber then
		return [ SEC`Zetas[2*k-1]*y : k in [1..SEC`Degree[1]] ];
	else
		return [y];
	end if;
end intrinsic;
					

function RS_AJMIntegrate( SEC, Ind, P )
// Integrate (naively) from P_i to P // This can be improved
	C<i> := SEC`ComplexField; C_0 := Zero(C); g := SEC`Genus;
	m := SEC`Degree[1]; n := SEC`Degree[2];
	DFF := &cat[ DFF_i : DFF_i in SEC`HolomorphicDifferentials ];

	a := SEC`BranchPoints[Ind];
	P_x := P[1];
	P_y := P[2];

	// Make CCV vector
	CCV := RS_MakeCCVectors(Ind,P_x,SEC`BranchPoints);
	UP := #[ z : z in CCV | Re(z) gt 0 ] mod 2;
	assert #CCV eq n-1;

	// Factors due to change of variable
	Fact1 := (P_x-a)/2;
	Fact2 := (P_x+a)/(P_x-a);
	
	vprint SE,2 : "################### Next Integral ###################";
	Integral := [ C_0 : j in [1..g] ];
	
	// Initiate on x = 0, dx = 1
	z := RS_NRootAC(0,CCV,SEC`Zetas,m);
	for j in [1..g] do
		w := DFF[j];
		Denom := z^w[2];
		Integral[j] +:= (Fact2^w[1]/Denom);
	end for;

	for t in [1..#SEC`Abscissas] do
		x := SEC`Abscissas[t];

		// Analytic continuation
		z1 := RS_NRootAC(x,CCV,SEC`Zetas,m);
		z2 := RS_NRootAC(-x,CCV,SEC`Zetas,m);
		
		// End point not a branch point
		c1 := (1-x)^(1/m);
		c2 := (1+x)^(1/m);

		Enum1 := (x + Fact2);
		Enum2 := (-x + Fact2);

		for j in [1..g] do
			w := DFF[j];
			Denom1 := z1^w[2];
			Denom2 := z2^w[2];
			dx := SEC`Weights[t][1] * SEC`Weights[t][2]^(2-(2*w[2]/m));
			Integral[j] +:= ( Enum1^w[1] * c1^w[2] / Denom1 + Enum2^w[1] * c2^w[2] / Denom2) * dx;
		end for;
	end for;

	// Correct sheet?
	P_y_AC := Exp( (n/m) * Log(P_x-a))  * SEC`Zetas[ UP + 1]  * RS_NRootAC(1,CCV,SEC`Zetas,m) / SEC`LeadingCoeff;

	// Shifting number
	s_ := (m/(2*SEC`Pi)) * ( Arg(P_y) - Arg(P_y_AC) );
	print "s:",s_;
	s := Round(s_);
	// Check: k \in \Z ?
	if Abs(s-s_) gt SEC`Error then
		print "CCV:",CCV;
		print "SEC`LeadingCoeff:",SEC`LeadingCoeff;
		print "RS_NRootAC(1,CCV,SEC`Zetas,m):",RS_NRootAC(1,CCV,SEC`Zetas,m);
		print "P_y_AC:",P_y_AC;
		print "Arg(P_y):",Arg(P_y);
		print "Arg(P_y_AC):",Arg(P_y_AC);
		print "Poly:",SEC`DefiningPolynomial;
		assert Abs(s-s_) lt 10^-10;
	end if;

	for j in [1..g] do
		w := DFF[j];
		Factor := Fact1^((w[1]+1)-(n*w[2]/m));
		Pow := -(UP +2*s)*w[2] mod (2*m) + 1;
		Integral[j] *:= SEC`LeadingCoeff^w[2] * SEC`Zetas[Pow] * SEC`StepLength * Factor;
	end for;
        return Matrix(C,g,1,Integral);
end function;


function RS_UpdateIntegrationParameters(SEC,AJMTau)
	if AJMTau lt SEC`Tau then
		print "OldTau:",SEC`Tau;
		print "NewTau:",AJMTau;
		SEC`Tau := AJMTau;
		SEC`Abscissas, SEC`Weights, SEC`StepLength := RS_SEDEIntegrationParameters(SEC`LowPrecBranchPoints,SEC`MST,AJMTau,SEC`Degree[1],SEC`ComplexField);
	end if;
	return 0;
end function;

intrinsic RS_InitiateAbelJacobi( SEC::SECurve : Recompute := false )
{ Compute the Abel-Jacobi map of the point P = (x,y) }

	if Recompute or not assigned SEC`AbelJacobi then

		// Complex field
		C<i> :=  SEC`ComplexField;
		// Real field
		R := SEC`RealField; Q := RationalField();

                // Compute period matrix
		RS_SEPM(SEC:Small:=true,Recompute:=Recompute);

		// Compute 'map' of the spanning tree
		RS_TreeMatrix(SEC:Recompute:=Recompute);

		// Compute reduction matrix
		RS_ReductionMatrix(SEC:Recompute:=Recompute);

		// Compute Abel-Jacobi map between ramification points
		RS_AJMRamPoints(SEC:Recompute:=Recompute);

		// Define Abel-Jacobi map
		AbelJacobi := function ( D )
		// Computes the Abel-Jacobi map of the divisor [D - deg(D)*P_0]
			assert Type(D) eq SEDivisor;
			assert D`DefiningPolynomial eq SEC`DefiningPolynomial;

			ComplexIntegral := Matrix(C,SEC`Genus,1,[]);
			RationalIntegral := Matrix(Q,2*SEC`Genus,1,[]);
			RealIntegral := Matrix(R,2*SEC`Genus,1,[]);
			
			AJMTau := SEC`Tau;
			ToBeComputed := [];
			for j in [1..#D`Support] do
				P := D`Support[j];
				v_P :=  D`Coefficients[j];
				if #P eq 3 then
					// Infinite points
					if SEC`Degree[3] eq 1 then
						RationalIntegral +:= v_P * SEC`AJM_InftyPoints[P[2]];
					else
						RealIntegral +:= v_P * SEC`AJM_InftyPoints[P[2]];
					end if;
				else
					// Finite points
					Dist, Ind := RS_Distance(P[1],SEC`BranchPoints);  
					RationalIntegral +:= v_P * SEC`AJM_RamPoints[Ind];
					if Dist gt SEC`Error then
						//print "Dist:",Dist; print "Ind:",Ind;
						AJMTau := Min(AJMTau,RS_AJMTau(C_20!P[1],Ind,SEC`LowPrecBranchPoints));
						Append(~ToBeComputed,<Ind,j>);
					end if;
				end if;
			end for;
			
			// Update integration parameters // rework this!
			ok := RS_UpdateIntegrationParameters(SEC,AJMTau);
		
			// Actual integrations from P_k to P
			for tp in ToBeComputed do
				ComplexIntegral +:= D`Coefficients[tp[2]] * RS_AJMIntegrate(SEC,tp[1],D`Support[tp[2]]);
			end for;

			// Special cases
			if ComplexIntegral eq Matrix(C,SEC`Genus,1,[]) then
				if RealIntegral eq Matrix(R,2*SEC`Genus,1,[]) then
					return Matrix(Q,2*SEC`Genus,1,[v - Round(v) : v in Eltseq(RationalIntegral) ]);
				else
					RealIntegral +:= ChangeRing(RationalIntegral,R);
					return Matrix(R,2*SEC`Genus,1,[v - Round(v) : v in Eltseq(RealIntegral) ]);
				end if;
			end if;

			// Reduce complex vector modulo (A,B)
			TotalRealIntegral := Matrix(R,2*SEC`Genus,1,RS_LatticeReduction(SEC,Eltseq(ComplexIntegral)));

			// Add integrals corresponding to ramification points and infinite points 
			TotalRealIntegral +:= ChangeRing(RationalIntegral,R) + RealIntegral;
			
			// Reduce again modulo \Z and return results
			//R_ := RealField(SEC`Prec);
			return Matrix(R,2*SEC`Genus,1,[v - Round(v) : v in Eltseq(TotalRealIntegral) ]);
		end function;

		// Save Abel-Jacobi map as attribute
		SEC`AbelJacobi := AbelJacobi;

		// Compute Abel-Jacobi map for points above infinity
		RS_AJMInftyPoints(SEC:Recompute:=Recompute);
	end if;	
end intrinsic;


intrinsic RS_TreeMatrix( SEC::SECurve : Recompute := false  )
{ Compute a matrix with paths in the spanning tree from P_0 -> P_i for all ramification points }
	if not assigned SEC`TreeMatrix or Recompute then
	// Compute period matrix
	RS_SEPM(SEC:Small:=false);
  	Points := SEC`BranchPoints;
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
	// Shift by basepoint for AJM
	P_0 := SEC`Basepoint; 
	P_0P_0 := TM[P_0];
	for j in [1..d] do
		TM[j] -:= P_0P_0;
	end for;
	SEC`TreeMatrix := TM;
	end if;
end intrinsic;


intrinsic RS_AJMRamPoints( SEC::SECurve : Recompute := false )
{ Abel-Jacobi map of ramification points }
	if not assigned SEC`AJM_RamPoints or Recompute then
	C<i> := ComplexField(SEC`Prec+10); Q := Rationals();
	g := SEC`Genus;
	m := SEC`Degree[1];
	n := SEC`Degree[2];
	RS_SEPM(SEC:Small:=false);
	RS_TreeMatrix(SEC);
	TM := SEC`TreeMatrix;
	AJMRamPoints := [];
	for j in [1..n] do
		V := Matrix(C,g,1,[]);
		TreePath := SEC`TreeMatrix[j];
		for k in [1..n-1] do
			V +:= TreePath[k] * SEC`ElementaryIntegrals[k];
		end for;
		V_red := RS_LatticeReduction(SEC,ElementToSequence(V));
		Q_red := [];
		for v in V_red do
			q := BestApproximation(v,m);
			assert Abs(q-v) lt SEC`Error;
			Append(~Q_red,q);
		end for;
		//Q_red := Matrix(Q,2*g,1,[ BestApproximation(v,m) : v in V_red ]);
		Append(~AJMRamPoints,Matrix(Q,2*g,1,Q_red));
	end for;
	SEC`AJM_RamPoints := AJMRamPoints;
	end if;
end intrinsic;
		

function RS_ComputeBranchPoints( SEC, C )
	f_x := ChangeRing(UnivariatePolynomial(Evaluate(SEC`DefiningPolynomial,[Parent(SEC`DefiningPolynomial).1,0])),C);
	return RS_Roots(f_x);
end function;


intrinsic RS_AJMInftyPoints( SEC::SECurve : Recompute:=false )
{ Abel-Jacobi map of points at infinity }
	if not assigned SEC`AJM_InftyPoints or Recompute then
	Q := Rationals();
	m := SEC`Degree[1];
	n := SEC`Degree[2];
	delta,a,b := Xgcd(m,n);
	SEC`AJM_InftyPoints := [];
	// Easy case
	if delta eq 1 then
		AJMInfty := a * &+[ V : V in SEC`AJM_RamPoints ];
		Append(~SEC`AJM_InftyPoints,Matrix(Q,2*SEC`Genus,1,[ z - Round(z) : z in Eltseq(AJMInfty)] ));
	else
	while a ge 0 do
		a -:= n;
		b +:= m;
	end while;
	assert a*m + b*n eq delta;
	M := Round(m/delta);
	N := Round(n/delta);

	// Update structure
	R := SEC`RealField;
	C<i> := ComplexField(SEC`Prec+Floor(SEC`Genus/2)+10);
	//C<i> := ComplexField(200);
	print "C:",C;
	Ct<t> := PolynomialRing(C); 
	SEC`Pi := Real(Pi(C));
	SEC`BranchPoints := RS_ComputeBranchPoints(SEC,C); 

	
	if M eq 1 then
		// Solve polynomial g(1,t) = 0 
		g := &*[ (1 - x*(t^M)) : x in SEC`BranchPoints ] - 1;
		assert Degree(g) in [n*M,(n-1)*M];
		g_C := Coefficients(g);
		Inf_Ord := Min([ j : j in [1..Degree(g)+1] | Abs(g_C[j]) gt SEC`Error ]) - 1;
		if Inf_Ord eq 1 then
			// Obtain the other integrals by multiplication with correct powers of zeta
			DFF := &cat[ DFF_i : DFF_i in SEC`HolomorphicDifferentials ];
			ZetaPows := [];
			for j in [1..SEC`Genus] do
				Pow := Round(-2*m*(a+b*N/M)*DFF[j][2]/delta) mod (2*m) + 1;
				Append(~ZetaPows,Pow);
			end for; 

			// Define g(1,t)/t^ord
			g_t := &+[ g_C[j+1] * t^(j-Inf_Ord) : j in [Inf_Ord..Degree(g)] | Abs(g_C[j+1]) gt SEC`Error ];
			// Compute roots
			Rts_g_t := Roots(g_t);
			// Get points and coefficients
			Points := < >; 
			Coeffs := [];
			for j in [1..#Rts_g_t] do
				t_j := Rts_g_t[j][1];
				x_j := 1/((t_j^M));
				y_j := 1/((t_j^N)*SEC`LeadingCoeff);
				Append(~Points,[x_j,y_j]);
				Append(~Coeffs,Rts_g_t[j][2]);
			end for;
			// Rational part
			RationalIntegral := Matrix(Q,2*SEC`Genus,1,[]);
			// Is zero a branch point?
			if Degree(g) eq (n-1)*M then
				Dist, Ind := RS_Distance(Zero(C),SEC`BranchPoints);
				assert Dist lt SEC`Error;
				RationalIntegral -:= (N*m-M) * SEC`AJM_RamPoints[Ind];
			end if;
			// Sort out ramification points // they occur here?!?
			ToBeComputed := []; AJMTau := SEC`Tau;
			for j in [1..#Points] do		
				Dist, Ind := RS_Distance(Points[j][1],SEC`BranchPoints);
				if Dist gt SEC`Error then
					RationalIntegral +:= Coeffs[j] * SEC`AJM_RamPoints[Ind];
					AJMTau := Min(AJMTau,RS_AJMTau(C_20!Points[j][1],Ind,SEC`LowPrecBranchPoints));
					Append(~ToBeComputed,<Ind,j>);
				else
					//error Error("Should not happen.");
					print "Dist:",Dist; print "Ind:",Ind;
					RationalIntegral +:= Coeffs[j] * SEC`AJM_RamPoints[Ind];
				end if;
			end for;

			// Update integration parameters // rework this!
			ok := RS_UpdateIntegrationParameters(SEC,AJMTau);

			// Actual integrations from P_k to P
			ComplexIntegrals := [ Matrix(SEC`ComplexField,SEC`Genus,1,[]) : j in [1..delta] ];
			for tp in ToBeComputed do
				ComplexIntegral0 := RS_AJMIntegrate(SEC,tp[1],Points[tp[2]]);
				ComplexIntegrals[1] +:= Coeffs[tp[2]] * ComplexIntegral0;
				for k in [2..delta] do
					CI0seq := Eltseq(ComplexIntegral0);
					ComplexIntegral0 := Matrix(SEC`ComplexField,SEC`Genus,1,[ SEC`Zetas[ZetaPows[j]] * CI0seq[j] : j in [1..SEC`Genus]]);
					ComplexIntegrals[k] +:= Coeffs[tp[2]] * ComplexIntegral0;
				end for;
			end for;

			// Reduce complex vector modulo (A,B)
			RealIntegrals := [ Matrix(R,2*SEC`Genus,1,RS_LatticeReduction(SEC,Eltseq(CI))) : CI in ComplexIntegrals ];
			// Add integrals corresponding to ramification points
			TotalRealIntegrals := [ ChangeRing(RationalIntegral,R) + RealIntegrals[j] : j in [1..delta] ];
			// Reduce again modulo \Z and return results
			SEC`AJM_InftyPoints := [ Matrix(R,2*SEC`Genus,1,[ Round(v)-v : v in Eltseq(TRI) ]) : TRI in TotalRealIntegrals ];
			//print "SEC`AJM_InftyPoints :",SEC`AJM_InftyPoints;
		end if;
	end if;
	if #SEC`AJM_InftyPoints eq 0 then
		for k in [0..delta-1] do
			r := Exp(2*SEC`Pi*i*k/delta);
			g_t := &*[ (1 - x*(r+t)^b*(t^M)) : x in SEC`BranchPoints ] - (r+t)^delta;
			g_C := Coefficients(g_t);
			assert Degree(g_t) in [n*(M+b),(n-1)*(M+b)];
			assert Abs(g_C[2]) gt SEC`Error;
			Rts_g_t := RS_Roots(g_t);
			assert #Rts_g_t eq Degree(g_t);
			Points := < >; Coeffs := [];
			for j in [1..Degree(g_t)] do
				t_j := Rts_g_t[j];
				if Abs(t_j) gt SEC`Error then
					x_j := 1/(r+t_j)^b/(t_j^M);
					y_j := (r+t_j)^a/((t_j^N)*SEC`LeadingCoeff);
					Append(~Points,[x_j,y_j]);
					Append(~Coeffs,1);
				end if;
			end for;
			assert &+Coeffs in [n*b+n*M-1,n*b+n*M-b-M-1];
			if Degree(g_t) eq (n-1)*(M+b) then
				Append(~Points,[0,0]);
				Append(~Coeffs,-N*m+M+b);
			end if;
			for j in [1..n] do
				Append(~Points,[SEC`BranchPoints[j],0]);
				Append(~Coeffs,-b);
			end for;
			D := RS_SEDivisor(Points,Coeffs,SEC:Check:=true);
			V := -SEC`AbelJacobi(D);
			Append(~SEC`AJM_InftyPoints,Matrix(R,2*SEC`Genus,1,[ z - Round(z) : z in Eltseq(V) ]));
		end for;
	end if;
	end if;
	Z := (m/delta) * &+[ V : V in SEC`AJM_InftyPoints ]; print "Z:",Z;
	assert &and[ Abs(z-Round(z)) lt SEC`Error : z in Eltseq(Z) ];
	end if;
end intrinsic;


intrinsic RS_SEPM( SEC::SECurve : Recompute := false, Small := true )
{ Computes period matrices associated to the superelliptic curve defined by f to precision of SEC }

	if Recompute or &and[Small,not assigned SEC`SmallPeriodMatrix] or &and[not Small, not assigned SEC`BigPeriodMatrix] then 

	// Degrees
	m := SEC`Degree[1]; n := SEC`Degree[2];

	// Branch points
	Points := SEC`BranchPoints;
	
	// Low precision branch points
	LowPrecPoints := [ C_20!Pt : Pt in Points ];

	// Complex fields and constants
	C_<i> := ComplexField(SEC`Prec+10);
	C<i> := SEC`ComplexField; C_0 := Zero(C);

	// Maximal spanning tree
	vprint SE,1 : "Constructing spanning tree...";
	if n ge 15 then
		vprint SE,1 : "using euclidean weights...";
		SEC`MST, SEC`Tau := RS_SEMST2(LowPrecPoints:ED:=true);
	else
		vprint SE,1 : "using holomorphic weights...";
		SEC`MST, SEC`Tau := RS_SEMST2(LowPrecPoints:ED:=false);
	end if;

	// Holomorphic differentials
	DFF := SEC`HolomorphicDifferentials;

	// Genus
	g := SEC`Genus;
	vprint SE,1 : "Genus:",g;

	// Integration parameters
	vprint SE,1 : "Computing Integration parameters...";
	SEC`Abscissas, SEC`Weights, SEC`StepLength := RS_SEDEIntegrationParameters(LowPrecPoints,SEC`MST,SEC`Tau,m,SEC`ComplexField);

	// Computing periods
	vprint SE,1 : "Integrating...";
	Integrals := [];
	SEC`ElementaryIntegrals := [];
	for k in [1..n-1] do
		vprint SE,2 : "Integrating edge nr.",k;
		I, EI := RS_SEIntegrate( SEC`MST[k],Points,DFF,SEC`Abscissas,SEC`Weights,SEC`StepLength,n,m,SEC`Zetas,SEC`LeadingCoeff );
		//I, EI := RS_SEIntegrate3( SEC`MST[k],Points,DFF,SEC`Abscissas,SEC`Weights,SEC`StepLength,n,m,Zetas,SEC`LeadingCoeff );
		Append(~SEC`ElementaryIntegrals,Matrix(C_,g,1,EI));
		Append(~Integrals,I);
	end for;

	// Period matrix
	PM_ := ZeroMatrix(C,g,(m-1)*(n-1));
	for k in [1..n-1] do
		for l in [1..m-1] do
			Ind := (m-1)*(k-1) + l;
			for j in [1..g] do
				PM_[j][Ind] := Integrals[k][j][l];
			end for;
		end for;
	end for;
	vprint SE,3: "Integrals:",PM_;
	SEC`PeriodMatrix := ChangeRing(PM_,C_);

	// Intersection matrix
	vprint SE,1 : "Computing spsm-matrix...";
	spsm_Matrix := [ [] : j in [1..n-1] ];
	for j in [1..n-1] do
		spsm_Matrix[j][j] := <1,m-1>;
		for l in [j+1..n-1] do
			spsm := RS_SEIntersection(SEC`MST[j],SEC`MST[l],LowPrecPoints,m,SEC`Zetas);
			vprint SE,2: "<sp,sm> pair of edges",SEC`MST[j],"and",SEC`MST[l],":",<spsm[1] mod m,spsm[2] mod m>;
			spsm_Matrix[j][l] := <spsm[1] mod m,spsm[2] mod m>;
		end for;
	end for;
	SEC`spsm_Matrix := spsm_Matrix;
	vprint SE,3: "spsm_Matrix:",spsm_Matrix;
	

	vprint SE,1 : "Computing intersection matrix...";
	// Building block matrices
	Blocks := [];
	// Block matrix for self-shifts build the diagonal of intersection matrix
	SelfShiftBlock := ZeroMatrix(Integers(),m-1,m-1);
	for Ind1 in [1..m-1] do
		for Ind2 in [Ind1+1..m-1] do
			if Ind1 + 1 - Ind2 mod m eq 0 then
				SelfShiftBlock[Ind1][Ind2] := 1;
				SelfShiftBlock[Ind2][Ind1] := -1;
			end if;
		end for;
	end for;

	// Build blocks for intersection matrix
	for j in [1..n-1] do
		Blocks[(j-1)*n+1] := SelfShiftBlock;
		for l in [j+1..n-1] do
			Block := ZeroMatrix(Integers(),m-1,m-1);
			if spsm_Matrix[j][l] ne <0,0> then
				for Ind1 in [1..m-1] do
					sp := spsm_Matrix[j][l][1];
					sm := spsm_Matrix[j][l][2];
					Ind2 := (Ind1 + sp) mod m;
					if Ind2 ne 0 then
						Block[Ind1][Ind2] := 1;
					end if;
					Ind2 := (Ind1 + sm) mod m;
					if Ind2 ne 0 then
						Block[Ind1][Ind2] := -1;
					end if;
				end for;
			end if;
			Blocks[(j-1)*(n-1)+l] := Block;
			Blocks[(l-1)*(n-1)+j] := -Transpose(Block);
		end for;
	end for;
	SEC`IntersectionMatrix := BlockMatrix(n-1,n-1,Blocks);
	assert Rank(SEC`IntersectionMatrix) eq 2*g;

	
	// Symplectic reduction of intersection matrix
	vprint SE,1 : "Performing symplectic reduction...";
	CF, ST := RS_SymplecticBasis(SEC`IntersectionMatrix);
	ST_C := ChangeRing(Transpose(ST),C);
	vprint SE,3: "ST:",ST;
	vprint SE,3: "CF:",CF;
	
	vprint SE,1 : "Matrix multiplication 1...";
	PMAPMB := PM_ * ST_C;
	
	PM_A := ColumnSubmatrixRange(PMAPMB,1,g);
	PM_B := ColumnSubmatrixRange(PMAPMB,g+1,2*g);
	vprint SE,3 : "Dependent columns:",ColumnSubmatrixRange(PMAPMB,2*g+1,Nrows(ST));
	vprint SE,3 : "PM_A:",PM_A;
	vprint SE,3 : "PM_B:",PM_B;

	// Compute big period matrix
	PM := HorizontalJoin(PM_B,PM_A) + ZeroMatrix(C_,g,2*g);
	SEC`BigPeriodMatrix := PM;

	if Small then
		vprint SE,1 : "Matrix inversion...";
		PM_AInv := PM_A^(-1);
		vprint SE,1 : "Matrix multiplication 2...";
		PM := PM_AInv * PM_B;
	
		// Smoothing out entries
		vprint SE,1 : "Smoothing out entries...";
		for j in [1..g] do
			for k in [1..g] do
				if Abs(Re(PM[j][k])) lt SEC`Error then
					PM[j][k] := Zero(C_) + Im(PM[j][k])*i;
				end if;
				if Abs(Im(PM[j][k])) lt SEC`Error then
					PM[j][k] := Re(PM[j][k]) + Zero(C_)*i;		
				end if;
			end for;
		end for;
		PM := ChangeRing(PM,C_);
		vprint SE,3 : "PM:",PM;


		// Testing for symmetry of the period matrix
		vprint SE,1 : "Testing symmetry...";
		MaxSymDiff := 0;
		for j in [1..g] do
			for k in [j+1..g] do
				MaxSymDiff := Max(MaxSymDiff,Abs(PM[j][k] - PM[k][j]));
			end for;
		end for;
		vprint SE,1 : "Maximal symmetry deviation:",MaxSymDiff;
	
		if MaxSymDiff ge SEC`Error then
			error Error("Period matrix not symmetric: Computation failed.");
		end if;	
		
		// Testing positive definiteness of the imaginary part of the period matrix
		vprint SE,1 : "Testing positive definite...";
		PM_Im := ZeroMatrix(RealField(SEC`Prec),g,g);
		for j in [1..g] do
			for k in [j..g] do
				PM_Im[j][k] := Real(Im(PM[j][k]));
				PM_Im[k][j] := PM_Im[j][k];
			end for;
		end for;
		assert IsPositiveDefinite(PM_Im);
		SEC`SmallPeriodMatrix := PM;
	end if;

	end if;
end intrinsic;


// Some functionality
intrinsic RS_DeleteAll( SEC::SECurve )
{ Deletes all computed data }

	delete SEC`MST;
	delete SEC`Tau;
	delete SEC`SmallPeriodMatrix;
	delete SEC`BigPeriodMatrix;
	delete SEC`AbelJacobi;
	delete SEC`Theta;
	delete SEC`spsm_Matrix;
	delete SEC`IntersectionMatrix;
	delete SEC`ElementaryIntegrals;

end intrinsic;


// Divisors
declare type SEDivisor;

declare attributes SEDivisor: DefiningPolynomial, Support, Coefficients, Degree;

// Constructor
intrinsic RS_SEDivisor(Points::Tup,Coeffs::SeqEnum[RngIntElt],SEC::SECurve : Check:=true) -> SEDivisor
{ Construct the divisor \sum v_P P }
	require &and[#Points eq #Coeffs,0 notin Coeffs] : "Every point in the support needs a non-zero coefficient.";
	SED := New(SEDivisor);		
	SED`Degree := &+Coeffs;
	SED`Coefficients := Coeffs;	
	SED`DefiningPolynomial := SEC`DefiningPolynomial;

	// Complex field
	C<i> :=  SEC`ComplexField; C_1 := One(C); C_0 := Zero(C);

	// Support
	Supp := < >;
	for P in Points do
		if #P eq 1 then
			x_P := C!P[1]; y_P := RS_PrincipalBranch(SEC,x_P)[1]; 
			Append(~Supp,[x_P,y_P]);
		elif #P eq 2 then
			x_P := C!P[1]; y_P := C!P[2];
			if Check then
				t := Evaluate(SEC`DefiningPolynomial,[x_P,y_P]);
				if Abs(t) gt 10^-10 then
					print "t:",t;
					print "x_P:",x_P; 
					print "y_P:",y_P;
					error Error("Not a point on the curve.");
				end if;
				
			end if;
			Append(~Supp,[C!P[1],C!P[2]]);
		elif #P eq 3 then
			z_P := C!P[3];
			if z_P eq C_0 then
				require &and[ P[1] eq 0,P[2] in [1..SEC`Degree[3]] ] :
					"Please enter points above infinity as [0,j,0], indexed by j in",[1..SEC`Degree[3]];
				Append(~Supp,[0,P[2],0]);
			else
				x_P := P[1]/z_P;
				y_P := P[2]/z_P;
				if Check then
					t := Evaluate(SEC`DefiningPolynomial,[x_P,y_P]);
					if Abs(t) gt 10^-10 then
						print "t:",t;
						print "x_P:",x_P; 
						print "y_P:",y_P;
						error Error("Not a point on the curve.");
					end if;
				end if;
				Append(~Supp,[C!P[1],C!P[2]]);
			end if;
		else
			require false : "Correct input: Points = < Type1,Type2,Type3,... >, where Type1 = [x_P] (just x-coordinate), Type2 = [x_P,y_P] (point on the curve), Type3 = [0,j,0] (point above infinity)"; 
		end if;
	end for;
	SED`Support := Supp;
	return SED;
end intrinsic;

// Printing
intrinsic Print( SED::SEDivisor )
{ Print Riemann surface }
	print " Divisor on curve given by 0 =",SED`DefiningPolynomial;
	print " Degree:",SED`Degree;
	print " Support:",SED`Support;
	print " Coefficients:",SED`Coefficients; 
end intrinsic;

intrinsic RS_SEZeroDivisor(Points::Tup,Coeffs::SeqEnum[RngIntElt],SEC::SECurve) -> SEDivisor
{ Creates a zero divisor on SEC by adding -Deg * P_0 }
	Deg := &+Coeffs;
	Append(~Points,[SEC`BranchPoints[SEC`Basepoint],0]);
	Append(~Coeffs,-Deg);
	return RS_SEDivisor(Points,Coeffs,SEC);
end intrinsic;	



// Not needed
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
function taured(nz,tau);
  g := NumberOfRows(tau);
  C := BaseRing(tau);
  Itau := Matrix(g,g,[Im(a) : a in ElementToSequence(tau)]);
  dum := Itau^-1*Matrix(g,1,[Im(zi) : zi in ElementToSequence(nz)]);
  v1 := Matrix(C,g,1,[Round(di) : di in ElementToSequence(dum)]);
  nz := nz - tau*v1;
  return nz - Matrix(C,g,1,[Round(Re(di)) : di in ElementToSequence(nz)]);
end function;





// Integrating to infty // not working properly // use moving lemma instead

function RS_NRootACInfty(x,p,Points,Zeta,N:long:=false)
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

intrinsic RS_AJInfinity( SEC::SECurve : Ind := 1 ) -> SeqEnum[FldComElt]
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

	InfTau := RS_SEInfTau( x_P, Points : Lambda := C_Pi/2 );
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
	z := RS_NRootACInfty(0,x_P,Points,SEC`Zetas,m);
	for j in [1..g] do
		w := DFF[j];
		Denom := z^w[2];
		Integral[j] +:= (1/Denom);
	end for;

	for t in [1..#SEC`Abscissas] do
		x := SEC`Abscissas[t];

		// Analytic continuation
		z1 := RS_NRootACInfty(x,x_P,Points,SEC`Zetas,m);
		z2 := RS_NRootACInfty(-x,x_P,Points,SEC`Zetas,m);
		
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


	P_y_AC := RS_NRootACInfty(-1,x_P,Points,SEC`Zetas,m); // * 2^(d/N)
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
end intrinsic;









