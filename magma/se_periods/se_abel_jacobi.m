/******************************************************************************

 Copyright (C) 2017 Christian Neurohr

 ******************************************************************************/

// Import global settings
//import "superelliptic.m": RS_SEIntegrate, RS_SEIntegrate3, AC_mthRoot, RS_SEMST, RS_SEMST2, RS_MakeCCVectors, RS_SETau3, RS_SEInfTau;
import "se_period_matrix.m": SE_PeriodMatrix;
import "se_de_int.m": DE_Int_Params_AbelJacobi;
import "se_help_funcs.m": MakeCCVector, Distance;
import "se_anal_cont.m": AC_mthRoot;

C20<i> := ComplexField(20);

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

intrinsic LatticeReduction( SEC::SECurve, V::SeqEnum[FldComElt] ) -> Mtrx
{ Wrapper }
	assert #V eq SEC`Genus;
	return LatticeReduction(SEC,[ Re(v) : v in V ] cat [ Im(v) : v in V ]);
end intrinsic;
intrinsic LatticeReduction( SEC::SECurve, V::SeqEnum[FldReElt] ) -> SeqEnum[FldReElt]
{ Reduce V \in \C^g modulo period matrix (A B) }
	R := SEC`RealField; g := SEC`Genus;
	// Assume V in \R^2g
	W := SEC`ReductionMatrix * Matrix(R,2*g,1,V);
	return [w - Round(w) : w in ElementToSequence(W)];
end intrinsic;

procedure SE_ReductionMatrix(SEC)
// Computes a matrix for reduction modulo big period matrix (A,B)
	if not assigned SEC`ReductionMatrix then
		if not assigned SEC`BigPeriodMatrix then
			SE_PeriodMatrix(SEC:Small:=false);
		end if;
		PM_AB := SEC`BigPeriodMatrix;
		print "PM_AB:",PM_AB;
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

intrinsic SE_PrincipalBranch( SEC::SECurve, x::FldComElt : Fiber := false ) -> SeqEnum[FldComElt]
{ Returns y(x) = f(x)^(1/m) }
	m := SEC`Degree[1];
	fx := SEC`LeadingCoeff * &*[ (x-P) : P in SEC`BranchPoints ];
	if m eq 2 then
		y := Sqrt(fx);
	else
		y := Exp( (1/m) * Log(fx) );
	end if;
	if Fiber then
		return [ SEC`Zetas[2*k-1]*y : k in [1..SEC`Degree[1]] ];
	else
		return [y];
	end if;
end intrinsic;
					

function SE_AJMIntegrate( SEC, Ind, P )
// Integrate (naively) from P_i to P // This can be improved
	C<i> := SEC`ComplexField; C_0 := Zero(C); g := SEC`Genus;
	m := SEC`Degree[1]; n := SEC`Degree[2];
	DFF := &cat[ DFF_i : DFF_i in SEC`HolomorphicDifferentials ];

	a := SEC`BranchPoints[Ind];
	P_x := P[1];
	P_y := P[2];

	// Make vector of centers of circumcircles
	CCV, up := MakeCCVector(Ind,P_x,SEC`BranchPoints);
	assert #CCV eq n+1;
	assert CCV[n] eq (P_x-a)/2;
	assert CCV[n+1] eq (P_x+a)/(P_x-a);
	assert &and[ Re(CCV[j]) gt 0 : j in [1..up] ];
	
	Integral := [ C_0 : j in [1..g] ];
	
	// Initiate on x = 0, dx = 1
	z := AC_mthRoot(0,CCV,SEC`Zetas,up,n-1,m);
	for j in [1..g] do
		w := DFF[j];
		Denom := z^w[2];
		Integral[j] +:= (CCV[n+1]^w[1]/Denom);
	end for;

	for t in [1..#SEC`Abscissas] do
		x := SEC`Abscissas[t];

		// Analytic continuation
		z1 := AC_mthRoot(x,CCV,SEC`Zetas,up,n-1,m);
		z2 := AC_mthRoot(-x,CCV,SEC`Zetas,up,n-1,m);
		
		// End point not a branch point
		c1 := (1-x)^(1/m);
		c2 := (1+x)^(1/m);

		Enum1 := (x + CCV[n+1]);
		Enum2 := (-x + CCV[n+1]);

		for j in [1..g] do
			w := DFF[j];
			Denom1 := z1^w[2];
			Denom2 := z2^w[2];
			dx := SEC`Weights[t][1] * SEC`Weights[t][2]^(2-(2*w[2]/m));
			Integral[j] +:= ( Enum1^w[1] * c1^w[2] / Denom1 + Enum2^w[1] * c2^w[2] / Denom2) * dx;
		end for;
	end for;

	// Correct sheet?
	P_y_AC := Exp( (n/m) * Log(P_x-a))  * SEC`Zetas[ up mod (2*m) + 1]  * AC_mthRoot(1,CCV,SEC`Zetas,up,n-1,m) / SEC`LeadingCoeff;

	// Shifting number
	s_ := (m/(2*SEC`Pi)) * ( Arg(P_y) - Arg(P_y_AC) );
	//print "s:",s_;
	s := Round(s_);
	// Check: k \in \Z ?
	if Abs(s-s_) gt SEC`Error then
		print "CCV:",CCV;
		print "SEC`LeadingCoeff:",SEC`LeadingCoeff;
		print "AC_mthRoot(1,CCV,SEC`Zetas,m):",AC_mthRoot(1,CCV,SEC`Zetas,up,n-1,m);
		print "P_y_AC:",P_y_AC;
		print "Arg(P_y):",Arg(P_y);
		print "Arg(P_y_AC):",Arg(P_y_AC);
		print "Poly:",SEC`DefiningPolynomial;
		assert Abs(s-s_) lt 10^-10;
	end if;

	for j in [1..g] do
		w := DFF[j];
		Factor := CCV[n]^((w[1]+1)-(n*w[2]/m));
		Pow := -(up +2*s)*w[2] mod (2*m) + 1;
		Integral[j] *:= SEC`LeadingCoeff^w[2] * SEC`Zetas[Pow] * SEC`StepLength * Factor;
	end for;
	/*print "Ind:",Ind;
	print "P:",P;
	print "Integral:",Integral;*/
        return Matrix(C,g,1,Integral);
end function;


intrinsic SE_InitiateAbelJacobi( SEC::SECurve : Recompute := false, Small := false )
{ Compute the Abel-Jacobi map of the point P = (x,y) }

	if Recompute or not assigned SEC`AbelJacobi then

		// Complex field
		C<i> :=  SEC`ComplexField;
		// Real field
		R := SEC`RealField; Q := RationalField();

                // Compute period matrix
		//SE_PeriodMatrix(SEC:Small:=Small,Recompute:=Recompute);

		// Compute 'map' of the spanning tree
		SE_TreeMatrix(SEC:Recompute:=Recompute);

		// Compute reduction matrix
		SE_ReductionMatrix(SEC);

		// Compute Abel-Jacobi map between ramification points
		SE_AJMRamPoints(SEC:Recompute:=Recompute);

		// Define Abel-Jacobi map
		AbelJacobi := function ( D )
		// Computes the Abel-Jacobi map of the divisor [D - deg(D)*P_0]
			assert Type(D) eq SEDivisor;
			assert D`DefiningPolynomial eq SEC`DefiningPolynomial;

			ComplexIntegral := Matrix(C,SEC`Genus,1,[]);
			RationalIntegral := Matrix(Q,2*SEC`Genus,1,[]);
			RealIntegral := Matrix(R,2*SEC`Genus,1,[]);
			
			ToBeComputed := [];
			Edges := [];
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
					Dist, Ind := Distance(P[1],SEC`BranchPoints);  
					RationalIntegral +:= v_P * SEC`AJM_RamPoints[Ind];
					if Dist gt SEC`Error then
						//print "Dist:",Dist; print "Ind:",Ind;
						Append(~ToBeComputed,<Ind,j>);
						Append(~Edges,<Ind,P[1]>);
					end if;
				end if;
			end for;
			
			// Compute integration parameters
			DEInt := DE_Int_Params_AbelJacobi( Edges, SEC );
		
			// Actual integrations from P_k to P
			for tp in ToBeComputed do
				ComplexIntegral +:= D`Coefficients[tp[2]] * SE_AJMIntegrate(SEC,tp[1],D`Support[tp[2]]);
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
			TotalRealIntegral := Matrix(R,2*SEC`Genus,1,LatticeReduction(SEC,Eltseq(ComplexIntegral)));

			// Add integrals corresponding to ramification points and infinite points 
			TotalRealIntegral +:= ChangeRing(RationalIntegral,R) + RealIntegral;
			
			// Reduce again modulo \Z and return results
			//R_ := RealField(SEC`Prec);
			return Matrix(R,2*SEC`Genus,1,[v - Round(v) : v in Eltseq(TotalRealIntegral) ]);
		end function;

		// Save Abel-Jacobi map as attribute
		SEC`AbelJacobi := AbelJacobi;

		// Compute Abel-Jacobi map for points above infinity
		SE_AJMInftyPoints(SEC:Recompute:=Recompute);
	end if;	
end intrinsic;


intrinsic SE_TreeMatrix( SEC::SECurve : Recompute := false  )
{ Compute a matrix with paths in the spanning tree from P_0 -> P_i for all ramification points }
	if not assigned SEC`TreeMatrix or Recompute then
	// Compute period matrix
	SE_PeriodMatrix(SEC:Small:=false);
  	Points := SEC`BranchPoints;
 	g := SEC`Genus; N := SEC`Degree[1]; d := SEC`Degree[2];
 
	TM := ZeroMatrix(Integers(),d,d-1);
	
	Taken := [ 0 : j in [1..d] ];
	Tree := SEC`SpanningTree`Edges;
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


intrinsic SE_AJMRamPoints( SEC::SECurve : Recompute := false )
{ Abel-Jacobi map of ramification points }
	if not assigned SEC`AJM_RamPoints or Recompute then
	//C<i> := ComplexField(SEC`Prec+10);
	C<i> := SEC`ComplexField;
	Q := Rationals();
	g := SEC`Genus;
	m := SEC`Degree[1];
	n := SEC`Degree[2];
	SE_PeriodMatrix(SEC:Small:=false);
	SE_TreeMatrix(SEC);
	TM := SEC`TreeMatrix;
	AJMRamPoints := [];
	for j in [1..n] do
		V := Matrix(C,g,1,[]);
		TreePath := SEC`TreeMatrix[j];
		for k in [1..n-1] do
			V +:= TreePath[k] * SEC`ElementaryIntegrals[k];
		end for;
		V_red := LatticeReduction(SEC,ElementToSequence(V));
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
		

function SE_ComputeBranchPoints( SEC, C )
	f_x := ChangeRing(SEC`ComplexPolynomial,C);
	Roots_f := Roots(f_x);
	return [ R[1] : R in Roots_f ];
end function;


intrinsic SE_AJMInftyPoints( SEC::SECurve : Recompute:=false )
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
	C<i> := ComplexField(SEC`Prec+2*SEC`Genus+10);
	//C<i> := ComplexField(50);
	//C<i> := SEC`ComplexField;
	print "C:",C;
	Ct<t> := PolynomialRing(C); 
	SEC`BranchPoints := SE_ComputeBranchPoints(SEC,C); 

	
	if M eq 1 then
		// Solve polynomial g(1,t) = 0 
		g := &*[ (1 - x*(t^M)) : x in SEC`BranchPoints ] - 1;
		assert Degree(g) in [n,n-1];
		g_C := Coefficients(g);
		Inf_Ord := Min([ j : j in [1..Degree(g)+1] | Abs(g_C[j]) gt SEC`Error ]) - 1;
		if Inf_Ord eq 1 then
			print "Good case!";
			// Obtain the other integrals by multiplication with correct powers of zeta
			jPows := SEC`HolomorphicDifferentials[4];
			ZetaPows := [];
			for k in [1..SEC`Genus] do
				Pow := Round(-2*m*(a+b*N)*jPows[k]/delta) mod (2*m) + 1;
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
				ltj := Log(t_j);
				x_j := Exp(-M*ltj);
				y_j := Exp(-N*ltj)/SEC`LeadingCoeff;
				//x_j := 1/((t_j^M));
				//y_j := 1/((t_j^N)*SEC`LeadingCoeff);
				Append(~Points,[x_j,y_j]);
				Append(~Coeffs,Rts_g_t[j][2]);
			end for;
			// Rational part
			RationalIntegral := Matrix(Q,2*SEC`Genus,1,[]);
			// Is zero a branch point?
			if Degree(g) eq (n-1)*M then
				Dist, Ind := Distance(Zero(C),SEC`BranchPoints);
				assert Dist lt SEC`Error;
				RationalIntegral -:= (N*m-M) * SEC`AJM_RamPoints[Ind];
			end if;
			// Sort out ramification points // they occur here?!?
			ToBeComputed := []; Edges := [];
			for j in [1..#Points] do		
				Dist, Ind := Distance(Points[j][1],SEC`BranchPoints);
				if Dist gt SEC`Error then
					RationalIntegral +:= Coeffs[j] * SEC`AJM_RamPoints[Ind];
					Append(~ToBeComputed,<Ind,j>);
					Append(~Edges,<Ind,Points[j][1]>);
				else
					//error Error("Should not happen.");
					print "Dist:",Dist; print "Ind:",Ind;
					RationalIntegral +:= Coeffs[j] * SEC`AJM_RamPoints[Ind];
				end if;
			end for;

			// Compute integration parameters
			DEInt := DE_Int_Params_AbelJacobi( Edges, SEC );
			

			// Actual integrations from P_k to P
			ComplexIntegrals := [ Matrix(SEC`ComplexField,SEC`Genus,1,[]) : j in [1..delta] ];
			for tp in ToBeComputed do
				ComplexIntegral0 := SE_AJMIntegrate(SEC,tp[1],Points[tp[2]]);
				ComplexIntegrals[1] +:= Coeffs[tp[2]] * ComplexIntegral0;
				for k in [2..delta] do
					CI0seq := Eltseq(ComplexIntegral0);
					ComplexIntegral0 := Matrix(SEC`ComplexField,SEC`Genus,1,[ SEC`Zetas[ZetaPows[j]] * CI0seq[j] : j in [1..SEC`Genus]]);
					ComplexIntegrals[k] +:= Coeffs[tp[2]] * ComplexIntegral0;
				end for;
			end for;

			// Reduce complex vector modulo (A,B)
			RealIntegrals := [ Matrix(R,2*SEC`Genus,1,LatticeReduction(SEC,Eltseq(CI))) : CI in ComplexIntegrals ];
			// Add integrals corresponding to ramification points
			TotalRealIntegrals := [ ChangeRing(RationalIntegral,R) + RealIntegrals[j] : j in [1..delta] ];
			// Reduce again modulo \Z and return results
			SEC`AJM_InftyPoints := [ Matrix(R,2*SEC`Genus,1,[ Round(v)-v : v in Eltseq(TRI) ]) : TRI in TotalRealIntegrals ];
			//print "SEC`AJM_InftyPoints :",SEC`AJM_InftyPoints;
		end if;
	end if;
	if #SEC`AJM_InftyPoints eq 0 then
		print "Bad case!";
		for k in [1..delta] do
			/*if M eq 1 and k eq delta then
				Append(~SEC`AJM_InftyPoints,-&+[ v : v in SEC`AJM_InftyPoints]);
				break;
			end if;*/
			r := Exp(2*Pi(C)*i*(k-1)/delta);
			//print "r:",r;
			g_t := &*[ (1 - x*(r+t)^b*(t^M)) : x in SEC`BranchPoints ] - (r+t)^delta;
			g_C := Coefficients(g_t);
			assert Degree(g_t) in [n*(M+b),(n-1)*(M+b)];
			assert Abs(g_C[2]) gt SEC`Error;
			Roots_gt := Roots(g_t);
			Rts_g_t := [ R[1] : R in Roots_gt ];
			assert #Rts_g_t eq Degree(g_t);
			Points := < >; Coeffs := [];
			for j in [1..Degree(g_t)] do
				t_j := Rts_g_t[j];
				if Abs(t_j) gt SEC`Error then
					lrtj := Log(r+t_j);
					ltj := Log(t_j);
					x_j := Exp(-((b*lrtj+M*ltj)));
					y_j := Exp(a*lrtj-N*ltj)/SEC`LeadingCoeff;
					//x_j := 1/(r+t_j)^b/(t_j^M);
					//y_j := (r+t_j)^a/t_j^N/SEC`LeadingCoeff;
					Append(~Points,[x_j,y_j]);
					Append(~Coeffs,1);
				end if;
			end for;
			//print "g_t:",g_t;
			assert &+Coeffs in [n*b+n*M-1,n*b+n*M-b-M-1];
			print "Have to compute:",&+Coeffs,"Integrals!";
			//print "Points:",Points;
			if Degree(g_t) eq (n-1)*(M+b) then
				Append(~Points,[0,0]);
				Append(~Coeffs,-N*m+M+b);
			end if;
			for j in [1..n] do
				Append(~Points,[SEC`BranchPoints[j],0]);
				Append(~Coeffs,-b);
			end for;
			D := SE_Divisor(Points,Coeffs,SEC:Check:=true);
			V := -SEC`AbelJacobi(D);
			Append(~SEC`AJM_InftyPoints,Matrix(R,2*SEC`Genus,1,[ z - Round(z) : z in Eltseq(V) ]));
		end for;
	end if;
	end if;
	Z := (m/delta) * &+[ V : V in SEC`AJM_InftyPoints ]; 
	for z in Eltseq(Z) do
		if Abs(z-Round(z)) gt SEC`Error then
			print "Z:",Z;
			print "SEC`DefiningPolynomial:",SEC`DefiningPolynomial;
			print "SEC`AJM_InftyPoints:",SEC`AJM_InftyPoints;
			assert false;
		end if;
	end for;
	//assert &and[ Abs(z-Round(z)) lt SEC`Error : z in Eltseq(Z) ];
	end if;
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
