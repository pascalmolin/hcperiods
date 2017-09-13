/*******************************************************************************

 Copyright (C) 2017 Christian Neurohr

*******************************************************************************/

// Import functions
import "se_de_int.m": DE_Integrals_Edge_AJM, DE_Int_Params;
import "se_abel_jacobi.m": PeriodLatticeReduction;
import "se_help_funcs.m": IsPoint,Distance;
import "se_spanning_tree.m": DE_Params_AJM;


intrinsic SE_AJM_InftyPoints( k::RngIntElt, SEC::SECurve )
{ Abel-Jacobi map of D = [ P_infty^k - P_0 ], where P_infty^k = (\zeta_\delta^k,0) in the (r,t)-model of the curve }
	m := SEC`Degree[1];
	n := SEC`Degree[2];
	delta,a,b := Xgcd(m,n);
	require k in [1..delta] : "There are only",delta,"infinite points.";
	Q := Rationals();
	R := SEC`RealField;
	if not IsDefined(SEC`AJM_InftyPoints,k) then 
		if delta eq 1 then
			Append(~SEC`AJM_InftyPoints,SEC`AJM_SumOfInftyPoints);
		elif &and[IsDefined(SEC`AJM_InftyPoints,j) : j in Remove([1..delta],k) ] then
			LastInftyPoint := SEC`AJM_SumOfInftyPoints - &+[ SEC`AJM_InftyPoints[j] : j in Remove([1..delta],k) ];
			SEC`AJM_InftyPoints[k] := Matrix(R,2*SEC`Genus,1,[ v-Round(v) : v in Eltseq(LastInftyPoint) ]);
		else
		while a ge 0 do
			a -:= n;
			b +:= m;
		end while;
		assert a*m + b*n eq delta;
		M := Round(m/delta);
		N := Round(n/delta);
		
		// Need more precision
		C<I> := ComplexField((Ceiling((m/10))+1)*SEC`Prec+2*n);
		Ct<t> := PolynomialRing(C); 

		// Embed data
		f_x := ChangeRing(SEC`DefiningPolynomial,C);
		Rts := Roots(f_x);
		BPs := [ R[1] : R in Rts ];

		// m-th root of leading coefficient
		LC_mi := Exp( (1/m) * Log(LeadingCoefficient(f_x)));
		if M eq 1 then
			// Solve polynomial g(1,t) = 0 
			g := &*[ (1 - x*(t^M)) : x in BPs ] - 1;
			assert Degree(g) in [n,n-1];
			g_C := Coefficients(g);
			Inf_Ord := Min([ j : j in [1..Degree(g)+1] | Abs(g_C[j]) gt SEC`Error ]) - 1;
			if Inf_Ord eq 1 then
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
					y_j := Exp(-N*ltj)*LC_mi;
					Append(~Points,<[x_j,y_j],Rts_g_t[j][2]>);
				end for;
				vprint SE,1 : "Good case: have to compute:",#Points,"Integrals!";
				// Rational part
				RationalIntegral := Matrix(Q,2*SEC`Genus,1,[]);
				// Is zero a branch point?
				if Degree(g) eq (n-1)*M then
					Dist, Ind := Distance(Zero(C),BPs);
					assert Dist lt SEC`Error;
					RationalIntegral -:= (N*m-M) * SEC`AJM_RamificationPoints[Ind];
				end if;
				// Sort out ramification points // they occur here?!?
				ComplexEdges := [];
				for P in Points	do
					Dist, Ind := Distance(P[1][1],BPs);
					RationalIntegral +:= P[2] * SEC`AJM_RamificationPoints[Ind];
					//if Dist gt SEC`Error then
					if Dist gt SEC`Error then
						Append(~ComplexEdges,Append(P,Ind));
					end if;
				end for;

				// Integration parameters
				Params, ComplexEdges := DE_Params_AJM(ComplexEdges,SEC);
				vprint SE,2 : "Parameter(AJM):",Params;
				ExtraPrec := 2*Ceiling(Log(10,Params[1]/SEC`SpanningTree`Params[1]));
				if ExtraPrec gt 0 then
					vprint SE,2 :"Extra precision (AJM):",ExtraPrec;
					C<I> := ComplexField(Precision(SEC`ComplexField)+ExtraPrec);
					SEC`ComplexField := C;
					// Update curve data
					f_x := ChangeRing(SEC`DefiningPolynomial,C);
					SEC`ComplexPolynomial := f_x;
					// Roots of f_x
					Roots_fx := Roots(f_x);
					SEC`BranchPoints := [ R[1] : R in Roots_fx ];
				end if;
				DEInts := DE_Integration(Params,SEC:AJM); NInts := #DEInts;

				// Actual integrations from P_k to P
				ComplexIntegrals := [ Matrix(SEC`ComplexField,SEC`Genus,1,[]) : j in [1..delta] ];
				for CE in ComplexEdges do
					l := NInts;
					while CE[4] lt DEInts[l]`r do
						l -:= 1;
					end while;
					ComplexIntegral0 :=  DE_Integrals_Edge_AJM(CE,SEC,DEInts[l]);
					ComplexIntegrals[1] +:= CE[2] * Matrix(SEC`ComplexField,SEC`Genus,1,ComplexIntegral0);
					for k in [2..delta] do
						CI0seq := Eltseq(ComplexIntegral0);
						ComplexIntegral0 := Matrix(SEC`ComplexField,SEC`Genus,1,[ SEC`Zetas[ZetaPows[j]] * CI0seq[j] : j in [1..SEC`Genus]]);
						ComplexIntegrals[k] +:= CE[2] * ComplexIntegral0;
					end for;
				end for;

				// Reduce complex vector modulo (A,B)
				RealIntegrals := [ Matrix(R,2*SEC`Genus,1,PeriodLatticeReduction(Eltseq(CI),SEC)) : CI in ComplexIntegrals ];
				// Add integrals corresponding to ramification points
				TotalRealIntegrals := [ ChangeRing(RationalIntegral,R) + RealIntegrals[j] : j in [1..delta] ];

				// Reduce again modulo \Z and return results in correct order (according to the paper)
				SEC`AJM_InftyPoints := [ Matrix(R,2*SEC`Genus,1,[ Round(v)-v : v in Eltseq(TotalRealIntegrals[k]) ]) : k in [2..delta] ];
				Append(~SEC`AJM_InftyPoints,Matrix(R,2*SEC`Genus,1,[ Round(v)-v : v in Eltseq(TotalRealIntegrals[1]) ]));
			end if;
		end if;
		if M gt 1 or Inf_Ord ne 1 then
			r := Exp(2*Pi(C)*I*k/delta);
			g_t := &*[ (1 - x*(r+t)^b*(t^M)) : x in BPs ] - (r+t)^delta;
			g_C := Coefficients(g_t);
			//g_t := &+[ g_C[j] * t^(j-1) : j in [1..Degree(g_t)+1] | Abs(g_C[j]) gt SEC`Error ];
			g_t := &+[ g_C[j] * t^(j-1) : j in [1..Degree(g_t)+1] ];
			assert Degree(g_t) in [n*(M+b),(n-1)*(M+b)];
			assert Abs(g_C[2]) gt SEC`Error;
			Roots_gt := Roots(g_t);
			Rts_g_t := [ R[1] : R in Roots_gt ];
			assert #Rts_g_t eq Degree(g_t);
			Points := [ ];
			for j in [1..Degree(g_t)] do
				t_j := Rts_g_t[j];
				if Abs(t_j) gt SEC`Error then
					lrtj := Log(r+t_j);
					ltj := Log(t_j);
					x_j := Exp(-((b*lrtj+M*ltj)));
					y_j := Exp(a*lrtj-N*ltj)*LC_mi;
					Append(~Points,<[x_j,y_j],1>);
				end if;
			end for;
			assert #Points in [n*b+n*M-1,n*b+n*M-b-M-1];
			vprint SE,1 : "Bad case: have to compute:",#Points,"Integrals!";
			if Degree(g_t) eq (n-1)*(M+b) then
				Append(~Points,<[0,0],N*m+M+b>);
			end if;
			// Compute divisor	
			D := SE_Divisor(Points,SEC:Check:=false);
			V := -SE_AbelJacobi(D,[BPs[SEC`Basepoint],0],SEC);
			// Substract the sum of infinite points
			V -:= b * ChangeRing(SEC`AJM_SumOfInftyPoints,R);
			SEC`AJM_InftyPoints[k] := Matrix(R,2*SEC`Genus,1,[ z - Round(z) : z in Eltseq(V) ]);
		end if;
		end if;
	end if;
end intrinsic;
