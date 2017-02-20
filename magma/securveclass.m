// Import global settings
import "superelliptic.m": RS_SEIntegrate, RS_SEIntegrate3, RS_NRootAC, RS_SEMST, RS_SEMST2, RS_MakeCCVectors, RS_SETau3, RS_SEInfTau;
import "comparefunctions.m": RS_CompareByFirstComplexEntry,RS_CompareFldComElt;
import "globalprecision.m": RS_Config;
import "pathmethods.m": RS_ChoiceOfPathStar;
import "cvector.m": RS_Quotrem_I, RS_Quotrem_II;

C_20<i> := ComplexField(20);
C_Pi := Real(Pi(C_20));

// Riemann surfaces as type RieSrf defined here 
declare type SECurve;

declare attributes SECurve: DefiningPolynomial, Genus, Degree, BranchPoints, MST,  Tau, ComplexField, RealField, HolomorphicDifferentials, Prec, Theta, SmallPeriodMatrix, BigPeriodMatrix, Pi, PeriodMatrix, ReductionMatrix, Abscissas, Weights, StepLength, AJM_RamPoints, AJM_Infinity, AbelJacobi, TreeMatrix, Error,
ElementaryIntegrals, spsm_Matrix, IntersectionMatrix, Basepoint, Zetas, LeadingCoeff;

// Constructor
intrinsic RS_SECurve( f::RngMPolElt : Prec := -1 ) -> SECurve
{ Creates the Riemann surface defined by f(x,y) = 0 and computes its function field and genus }
	
	// Requirements
	m := Degree(f,2); n := Degree(f,1);
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

	// Complex field
	CC<i> := ComplexField(SEC`Prec+10);
	SEC`ComplexField := CC; SEC`Pi := Real(Pi(CC));

	// Real field
	SEC`RealField := RealField(SEC`Prec+5);

	// Degree of cover
	SEC`Degree := [m,n,Gcd(m,n)];	

	// Branch points
	f_x := ChangeRing(UnivariatePolynomial(Evaluate(f,[Parent(f).1,0])),CC);
	SEC`BranchPoints :=  RS_Roots(f_x);

	// Leading coefficient
	SEC`LeadingCoeff := (-LeadingCoefficient(f_x))^(1/m);

	// Holomorphic differentials and genus
	SEC`HolomorphicDifferentials := RS_SEDifferentials(m,n:Int3:=false);

	// Genus
	SEC`Genus := Round((1/2) * ((m-1)*(n-1) - SEC`Degree[3] + 1));

	// Root of unity powers
	SEC`Zetas := [ RS_ReducePrecision(Exp(k*SEC`Pi*i/m),5) : k in [0..2*m-1] ];

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
	print " Basis of holomorphic differentials:", assigned SEC`HolomorphicDifferentials;
	print " Period matrix C:", assigned SEC`PeriodMatrix;
	print " Big period matrix (A,B):",assigned SEC`BigPeriodMatrix;
	print " Small period matrix A^(-1)B:",assigned SEC`SmallPeriodMatrix;
	print " Theta function:",assigned SEC`Theta;
	print " Abel-Jacobi map:", assigned SEC`AbelJacobi;
	print " AJM between ramification points:", assigned SEC`AJM_RamPoints;
	print " AJM to infinity:", assigned SEC`AJM_Infinity;
	print " Basepoint:", assigned SEC`Basepoint;
	print " Elementary integrals:", assigned SEC`ElementaryIntegrals;
	print " Tree matrix:", assigned SEC`TreeMatrix;
	print " Reduction matrix:", assigned SEC`ReductionMatrix;
end intrinsic;


function RS_ChoosePath( SEC, x )
// Choose path such that tau for integration is maximal 
	BestTau := 0; BestInd := 0;
	Tau := SEC`Tau; n := SEC`Degree[2]; Points := SEC`BranchPoints;
	for j in [1..n] do
		Tau_k := 4.0; // > Pi/2
		for k in [1..n] do
			if k ne j then
				Tau_k := Min(Tau_k,RS_SETau3(x,Points[j],Points[k]));
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
		SEC`Abscissas, SEC`Weights, SEC`StepLength := RS_SEDEIntegrationParameters(Points,SEC`MST,BestTau,SEC`Degree[1],SEC`ComplexField);
	end if;
	return BestInd;	
end function;

intrinsic RS_LatticeReduction( SEC::SECurve, V::SeqEnum[FldComElt] ) -> Mtrx
{ Wrapper }
	assert #V eq SEC`Genus;
	return RS_LatticeReduction(SEC,[ Re(v) : v in V ] cat [ Im(v) : v in V ]);
end intrinsic;
intrinsic RS_LatticeReduction( SEC::SECurve, V::SeqEnum[FldReElt] ) -> Mtrx
{ Reduce V \in \C^g modulo period matrix (A B) }
	R := SEC`RealField; g := SEC`Genus; 
	if not assigned SEC`ReductionMatrix then
		if not assigned SEC`BigPeriodMatrix then
			RS_SEPM(SEC:Small:=false);
		end if;
		PM_AB := SEC`BigPeriodMatrix;
		M := Matrix(R,2*g,2*g,[]);
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
	// Assume V in \R^2g
	W := SEC`ReductionMatrix * Matrix(R,2*g,1,V);
	Z := Matrix(R,2*g,1,[ w - Round(w) : w in ElementToSequence(W) ]);
	return Z;
end intrinsic;



intrinsic RS_PrincipalBranch( SEC::SECurve, x::FldComElt : Fiber := false ) -> SeqEnum[FldComElt]
{ Returns y(x) = f(x)^(1/m) }
	y := SEC`LeadingCoeff * (&*[ x-SEC`BranchPoints[k] : k in [1..SEC`Degree[2]] ])^(1/SEC`Degree[1]);
	if Fiber then
		return [ SEC`Zetas[2*k-1]*y : k in [1..SEC`Degree[1]] ];
	else
		return [y];
	end if;
end intrinsic;
					

function RS_AJMIntegrate( SEC, Ind, P )
// Integrate (naively) from P_i to P // This can be improved
	C<i> := SEC`ComplexField; C_0 := Zero(C); g := SEC`Genus;
	m := SEC`Degree[1]; n := SEC`Degree[2]; Points := SEC`BranchPoints; 
	DFF := &cat[ DFF_i : DFF_i in SEC`HolomorphicDifferentials ];

	a := Points[Ind];
	P_x := P[1];
	P_y := P[2];
	
	// Make CCV vector
	CCV := RS_MakeCCVectors(Ind,P_x,Points);
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
	P_y_AC := Exp( (n/m) * Log(P_x-a))  * SEC`Zetas[ UP + 1]  * RS_NRootAC(1,CCV,SEC`Zetas,m);
	k_ := (m/(2*SEC`Pi)) * ( Arg(P_y) - Arg(P_y_AC) );
	k := Round(k_);
	// Check: k \in \Z ?
	//assert Abs(k - k_) lt 10^-10;
	if Abs(k-k_) gt 10^-10 then
		print "k_:",k;
		print "f:",SEC`DefiningPolynomial;
	end if;
	for j in [1..g] do
		w := DFF[j];
		Factor := Fact1^((w[1]+1)-(n*w[2]/m));
		Pow := -(UP +2*k)*w[2] mod (2*m) + 1;
		Integral[j] *:= (1/SEC`LeadingCoeff)^w[2] * SEC`Zetas[Pow] * SEC`StepLength * Factor;
	end for;
        return Matrix(C,g,1,Integral);
end function;


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

intrinsic RS_AJInfinity( SEC::SECurve ) -> SeqEnum[FldComElt]
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
	
	P := [x_P,Y_P[3]];
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

	k_ := (m/(2*SEC`Pi)) * ( Arg(P[2]) - Arg(P_y_AC) );
	print "k_:",k_;
	k := Round(k_);
	assert Abs(k - k_) lt 10^-10; // Check: k \in \Z ?
	for j in [1..g] do
		w := DFF[j];
		Pow := -2*w[2]*k mod (2*m) + 1;
		print "Pow:",Pow;
		Integral[j] *:= (1/SEC`LeadingCoeff)^w[2] * SEC`Zetas[Pow] * SEC`StepLength * Fact^(w[1]+1);
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

	Res := [ ReducedIntegral[j][1] : j in [1..2*g] ];

	return Res;
end intrinsic;

intrinsic RS_InitiateAbelJacobi( SEC::SECurve : Recompute := false )
{ Compute the Abel-Jacobi map of the point P = (x,y) }

	if Recompute or not assigned SEC`AbelJacobi then

		// Complex field
		C<i> :=  SEC`ComplexField;

                // Compute period matrix
		RS_SEPM(SEC:Small:=true,Recompute:=Recompute);

		// Compute 'map' of the spanning tree
		RS_TreeMatrix(SEC:Recompute:=Recompute);

		// Compute Abel-Jacobi map between ramification points
		RS_AJMRamPoints(SEC:Recompute:=Recompute);

		// Define Abel-Jacobi map
		AbelJacobi := function ( P )
		// Computes the Abel-Jacobi map of the divisor D = [P - P_0]
			if #P eq 1 then
				x_P := C!P[1];
				y_P := RS_PrincipalBranch(SEC,x_P)[1];
			elif #P eq 2 then
				x_P := C!P[1];
				y_P := C!P[2];
			elif #P eq 3 then
				z_P := C!P[3];
				if z_P eq Zero(C) then
					if SEC`Degree[3] eq 1 then
						return SEC`AJM_Infinity;
					else
						error Error("TBD.");
					end if;
				else
					x_P := P[1]/z_P;
					y_P := P[2]/z_P;
				end if;
			else
				error Error("Invalid input.");
			end if;
			print "x_P,y_P:",x_P,y_P;
			t := Evaluate(SEC`DefiningPolynomial,[x_P,y_P]);
			if Abs(t) gt 10^-10 then
				print "t:",t;
				error Error("Not a point on the curve.");
			end if;

			// Is x_P a branch point?
			Dist, Ind := RS_Distance(x_P,SEC`BranchPoints);
			print "Dist:",Dist; print "Ind:",Ind;	
			if Dist lt SEC`Error then
				return SEC`AJM_RamPoints[Ind];
			end if;

			// Find closest branch points // Update integration parameters
			Ind := RS_ChoosePath(SEC, x_P);

			// Integrate from P_Ind to P
			Integrals := RS_AJMIntegrate(SEC,Ind,<x_P,y_P>);
						
			// Reduce modulo lattice
			RedIntegrals := RS_LatticeReduction(SEC,ElementToSequence(Integrals));

			// Add integral from P_0 to P_Ind
			RedIntegrals +:= ChangeRing(SEC`AJM_RamPoints[Ind],SEC`RealField);

			// Reduce again
			Result := Matrix(SEC`RealField,2*SEC`Genus,1,[v - Round(v) : v in ElementToSequence(RedIntegrals) ]);
			
			return Result;
		end function;
		// Save function;
		SEC`AbelJacobi := AbelJacobi;
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
	// Shift by real basepoint
	P_0 := SEC`Basepoint; 
	P_0P_0 := TM[P_0];
	for j in [1..d] do
		TM[j] -:= P_0P_0;
	end for;
	SEC`TreeMatrix := TM;
	end if;
end intrinsic;


intrinsic RS_AJMRamPoints( SEC::SECurve : Recompute := false )
{ Test for Abel-Jacobi map of superelliptic curves }
	if not assigned SEC`AJM_RamPoints or Recompute then
	C<i> := ComplexField(SEC`Prec); Q := Rationals();
	g := SEC`Genus;
	m := SEC`Degree[1];
	n := SEC`Degree[2];
	delta,mu,nu := Xgcd(m,n);
	RS_SEPM(SEC:Small:=false);
	RS_TreeMatrix(SEC);
	TM := SEC`TreeMatrix;
	AJMRamPoints := [];
	for j in [1..n] do
		V := Matrix(C,g,1,[]);
		TreePath := SEC`TreeMatrix[j];
		for j in [1..n-1] do
			V +:= TreePath[j] * SEC`ElementaryIntegrals[j];
		end for;
		V_red := RS_LatticeReduction(SEC,ElementToSequence(V));
		Q_red := Matrix(Q,2*g,1,[ BestApproximation(v,m) : v in ElementToSequence(V_red) ]);
		Append(~AJMRamPoints,Q_red);
	end for;
	SEC`AJM_RamPoints := AJMRamPoints;
	if delta eq 1 then
		AJMInfinity_nonred := mu * &+[ V : V in AJMRamPoints ];
		AJMInfinity := Matrix(Q,2*g,1,[ z - Round(z) : z in ElementToSequence(AJMInfinity_nonred) ]);
		SEC`AJM_Infinity := AJMInfinity;
		print "mu:",mu;
		print "AJMRamPoints:",AJMRamPoints;
		print "AJMInfinity:",AJMInfinity;
		print "m*AJMInfinity:",m*AJMInfinity;
	end if;
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
	C_<i> := ComplexField(SEC`Prec);
	C<i> := SEC`ComplexField; C_0 := Zero(C);

	// Maximal spanning tree
	vprint SE,1 : "Constructing spanning tree...";
	t := Cputime();
	if n ge 15 then
		vprint SE,1 : "using euclidean weights...";
		SEC`MST, SEC`Tau := RS_SEMST2(LowPrecPoints:ED:=true);
	else
		vprint SE,1 : "using holomorphic weights...";
		SEC`MST, SEC`Tau := RS_SEMST2(LowPrecPoints:ED:=false);
	end if;
	Cputime(t);

	// Holomorphic differentials
	DFF := SEC`HolomorphicDifferentials;

	// Genus
	g := SEC`Genus;
	vprint SE,1 : "Genus:",g;

	// Integration parameters
	vprint SE,1 : "Computing Integration parameters...";
	t := Cputime();
	SEC`Abscissas, SEC`Weights, SEC`StepLength := RS_SEDEIntegrationParameters(LowPrecPoints,SEC`MST,SEC`Tau,m,SEC`ComplexField);
	Cputime(t);

	// Computing periods
	vprint SE,1 : "Integrating...";
	t := Cputime();
	Integrals := [];
	SEC`ElementaryIntegrals := [];
	for k in [1..n-1] do
		vprint SE,2 : "Integrating edge nr.",k;
		I, EI := RS_SEIntegrate( SEC`MST[k],Points,DFF,SEC`Abscissas,SEC`Weights,SEC`StepLength,n,m,SEC`Zetas,SEC`LeadingCoeff );
		//I, EI := RS_SEIntegrate3( SEC`MST[k],Points,DFF,SEC`Abscissas,SEC`Weights,SEC`StepLength,n,m,Zetas,SEC`LeadingCoeff );
		Append(~SEC`ElementaryIntegrals,Matrix(C_,g,1,EI));
		Append(~Integrals,I);
	end for;
	Cputime(t);


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
	t := Cputime();
	spsm_Matrix := [ [] : j in [1..n-1] ];
	for j in [1..n-1] do
		spsm_Matrix[j][j] := <1,m-1>;
		for l in [j+1..n-1] do
			spsm := RS_SEIntersection(SEC`MST[j],SEC`MST[l],LowPrecPoints,m,SEC`Zetas);
			vprint SE,2: "<sp,sm> pair of edges",SEC`MST[j],"and",SEC`MST[l],":",<spsm[1] mod m,spsm[2] mod m>;
			spsm_Matrix[j][l] := <spsm[1] mod m,spsm[2] mod m>;
		end for;
	end for;
	Cputime(t);
	SEC`spsm_Matrix := spsm_Matrix;
	vprint SE,3: "spsm_Matrix:",spsm_Matrix;
	

	vprint SE,1 : "Computing intersection matrix...";
	t := Cputime();
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
	Cputime(t);
	assert Rank(SEC`IntersectionMatrix) eq 2*g;

	
	// Symplectic reduction of intersection matrix
	vprint SE,1 : "Performing symplectic reduction...";
	t := Cputime();
	CF, ST := RS_SymplecticBasis(SEC`IntersectionMatrix);
	ST_C := ChangeRing(Transpose(ST),C);
	Cputime(t);
	vprint SE,3: "ST:",ST;
	vprint SE,3: "CF:",CF;
	
	
	vprint SE,1 : "Matrix multiplication 1...";
	t := Cputime();
	PMAPMB := PM_ * ST_C;
	Cputime(t);
	
	PM_A := ColumnSubmatrixRange(PMAPMB,1,g);
	PM_B := ColumnSubmatrixRange(PMAPMB,g+1,2*g);
	vprint SE,3 : "Dependent columns:",ColumnSubmatrixRange(PMAPMB,2*g+1,Nrows(ST));
	vprint SE,3 : "PM_A:",PM_A;
	vprint SE,3 : "PM_B:",PM_B;

	
	// Compute big period matrix
	PM := HorizontalJoin(PM_B,PM_A) + ZeroMatrix(C_,g,2*g); // Dimensions???
	vprint SE,1 : "Period matrix:";
	SEC`BigPeriodMatrix := PM;


	if Small then
		vprint SE,1 : "Matrix inversion...";
		t := Cputime();
		PM_AInv := PM_A^-1;
		Cputime(t);
		vprint SE,1 : "Matrix multiplication 2...";
		t := Cputime();
		PM := PM_AInv * PM_B;
		Cputime(t);
	
		// Smoothing out entries
		vprint SE,1 : "Smoothing out entries...";
		t := Cputime();
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
	
		if MaxSymDiff ge SEC`Error then
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
	delete SEC`SmallPeriodMatrix;
	delete SEC`BigPeriodMatrix;
	delete SEC`AbelJacobi;
	delete SEC`Theta;
	delete SEC`spsm_Matrix;
	delete SEC`IntersectionMatrix;
	delete SEC`ElementaryIntegrals;

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

