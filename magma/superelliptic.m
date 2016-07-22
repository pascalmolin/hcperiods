//Import global settings
import "globalsettings.m": RS_Config;
import "complexgeo.m": RS_ComplexCompare;
import "integration.m": RS_IntParameters;
import "analyticcontinuation.m": RS_ACRecursion;
import "pathmethods.m":RS_FindSet;
import "comparefunctions.m":RS_CompareFldComElt;

declare verbose SE,3;
declare verbose SEP,2;

function RS_SETau3( P1,P2,P3 )
// Computes the maximal area of holomorphy between P1,P2 and P3
	Lambda := RS_GetGlobalPi()/2;
	CC := (2*P3 - P2 - P1) / (P2 - P1); // Circumcenter of the triangle defined by P1, P2, P3
	Tau_CC := Argsinh(Argtanh(CC)/Lambda);
	Tau_3 := Abs(Im(Tau_CC));
	return Tau_3;
end function;
procedure RS_SETau( Points, ~Edge )
// Computes Tau for the line segment Gamma
	for k in [1..#Points] do
		if &and[k ne Edge[2], k ne Edge[3]] then
			Edge[1] := Max(Edge[1],-RS_SETau3(Points[Edge[2]],Points[Edge[3]],Points[k]));
		end if;
	end for;
end procedure;
function RS_SEMST( Points )
// Computes a minimal spanning tree w.r.t holomorphicy
	//C_20 := ComplexField(20); 
	d := #Points;
	//LowPrecPts := [ C_20!Pt : Pt in Points ]; 
	Edges := [];
	
	for k in [1..d] do
		for l in [k+1..d] do
			if Re(Points[k]) le Re(Points[l]) then
				Edge := <-4.0,k,l>;
			else
				Edge := <-4.0,l,k>;
			end if;
			RS_SETau( Points, ~Edge );
			Append(~Edges,Edge);
		end for;
	end for;
	
	// Sort by first entry
	Sort(~Edges);
	Tau := -4.0;
	
	// Construct MST
	Sets := [ {v} : v in [1..d] ];
	MST_Edges := []; j := 1;
	
	while #MST_Edges lt d-1 do
		S1, s1 := RS_FindSet(Edges[j][2],Sets);
		S2, s2 := RS_FindSet(Edges[j][3],Sets);
		if S1 ne S2 then
			Append(~MST_Edges,<Edges[j][2],Edges[j][3]>);
			//print "Tau of edge nr.",j,":",Edges[j][1];
			Tau := Max(Tau,Edges[j][1]);
			Sets[s1] := Sets[s1] join Sets[s2];
			Remove(~Sets,s2);
		end if;
		j +:= 1;
	end while;
	assert #MST_Edges eq d-1;
	// Set Tau for integration
	//RS_SetIntTau((19/20)*(-Tau));
	//print "-Tau:",-Tau;
	return MST_Edges, -Tau;
end function;


function RS_ImSgn( z )
	if Im(z) eq 0 then
		return 1; // Need negative zero; example: g1
		//return 1;
	else
		return Sign(Im(z));
	end if;
end function;
function RS_NRootAC(x,PPP,Zeta,N:long:=false)
// Analytic continuation of the y = f(x)^(1/N) on one sheet
	CC<i> := RS_GetGlobalComplexField_Max();
	Pi := RS_GetGlobalPi();
	c := Exp(Pi*i/N);
	//print "c:",c;
	long := true;
	if #PPP eq 1 or long then
		prod := 1;
		for k in [1..#PPP] do
			prod *:= (x - PPP[k])^(1/N);
		end for;
		return prod;
	end if;
	z := x - PPP[1]; 
	Factor := 1;
	if Im(z) eq 0  and Re(z) lt 0 then
		z := -z;
		Factor *:= c;
	end if;
	isgnz := RS_ImSgn(z);
	for k in [2..#PPP] do
		z_ := x - PPP[k];
		if Im(z_) eq 0 and Re(z_) lt 0 then
			z_ := -z_;
			Factor *:= c;			
		end if;
		isgnz_ := RS_ImSgn(z_);
		z *:= z_;
		if isgnz eq isgnz_ then
			isgnz := RS_ImSgn(z);
			if isgnz ne isgnz_ then
				if isgnz_ gt 0 then
					Factor *:= Zeta;
				else
					Factor *:= (1/Zeta);
				end if;
			end if;
		else
			isgnz := RS_ImSgn(z);
		end if;
	end for;
	z := Factor * z^(1/N);
	/*
	zz := &*[ (x - PPP[k])^(1/N) : k in [1..#PPP] ];
	if Abs(z-zz) gt 10^-10 then
		print "x:",x;
		print "z:",z;
		print "zz:",zz;
		print "PPP:",PPP;
		print "N:",N;
		error Error("!!!!!!111");
	end if;
	*/
	return z;
end function;


function RS_RemovePoints(Points,List)
	// Remove Points of List
	assert #Points in [1,2];
	assert #Points ge 2;
	k := Position(List,Points[1]);
	l := Position(List,Points[2]);
	assert k ne l;
	print "l:",l;
	print "k:",k;
	if k eq 0 then
		PP := Remove(List,l);
	elif k lt l then
		PP := Remove(Remove(List,l),k);
	else
		PP := Remove(Remove(List,k),l);
	end if;
	return PP;
end function;


intrinsic RS_SEIntegration( f::RngMPolElt, Gamma::RSPath, Points::SeqEnum[FldComElt], DFF::SeqEnum[Tup] ) -> SeqEnum[FldComElt]
{ Integration on superelliptic curves }
	

	"##################### RS_SEIntegration ###################";
	// Variables
	d := Degree(f,1); assert d eq #Points;
	N := Degree(f,2); g := #DFF;
	assert #Points ge 2;

	a := Gamma`StartPt;
	b := Gamma`EndPt;

	PP := RS_RemovePoints([a,b],Points);
	PPP := [ (2*PP[k]-b-a)/(b-a) : k in [1..#PP] ]; 
	if #PP eq d-2 then
		cc := 0;
	else
		cc := 1;
	end if;

	Fact1 := (b-a)/2; print "Fact1:",Fact1;
	Fact2 := (b+a)/(b-a);

	
	RS_SetGlobalZeta(N);
	C<i> := RS_GetGlobalComplexField_Max(); Pi := RS_GetGlobalPi(); C_1 := One(C); CC_0 := Zero(C);
	Zeta := RS_GetGlobalZeta();

	print "Gamma:",Gamma;
	print "PP:",PP;
	
	RS_TauPath(PP,Gamma);
	RS_SetIntTau(Gamma`Tau);
	print "Gamma`Tau:",Gamma`Tau;
	
	Abscissas, Weights, StepLength := RS_SEIntegrationParameters(RS_GetIntTau(),N);
	
	//print "Abscissas:",Abscissas;
	//print "Weights:",Weights;
	//print "StepLength:",StepLength;

	// Factors due to change of variable
	Fact1 := (b-a)/2;
	Fact2 := (b+a)/(b-a);
	//IntFact := Exp(Pi*i/N) * StepLength;

	Integrals0 := [ CC_0 : j in [1..g] ];
	for t in [1..#Abscissas] do
		x := Abscissas[t];
		Enum := (x + Fact2);
		z := RS_NRootAC(x,PPP,Zeta,N);
		for j in [1..g] do
			w := DFF[j];
			Denom := z^w[2];
			dx := Weights[t][1] * Weights[t][2]^(2-(2*w[2]/N));
			Integrals0[j] +:= (Enum^w[1]/Denom) * dx;
		end for;
	end for;
	Integrals := [];
	for j in [1..g] do
		w := DFF[j];
		Factor := Fact1^((w[1]+1)-(d*w[2]/N));
		//Integrals0[j] *:= (1 - Zeta^(-w[2])) * IntFact * Factor;
		//Integrals0[j] *:= (1 - Zeta^(-w[2])) * Exp(Pi*i*w[2]/N) * StepLength * Factor;
		//Integrals0[j] *:= (1 - Zeta) * Exp(Pi*i/N) * StepLength * Factor;
		//Integrals0[j] *:= (1 - Zeta^(N-1)) * IntFact * Factor;

		Integrals0[j] *:= (Exp(-Pi*i*w[2]/N)-Exp(Pi*i*w[2]/N)) * StepLength * Factor;
		
		Append(~Integrals,[ Zeta^(l*w[2]) * Integrals0[j] : l in [0..N-1]]);
		//Append(~Integrals,[ -Zeta^l * Integrals0[j] : l in [0..N-1]]);
	end for;

	return Integrals;

end intrinsic;






function RS_SEIntegrate( Edge,Points,DFF,Abscissas,Weights,StepLength,d,N )

	CC<i> := RS_GetGlobalComplexField_Max(); Pi := RS_GetGlobalPi(); CC_0 := Zero(CC); CC_1 := One(CC); g := #DFF;
	Zeta := RS_GetGlobalZeta();
	a := Points[Edge[1]];
	b := Points[Edge[2]];
	assert Re(a) le Re(b);
	/*
	print "####################";
	print "For Edge:",Edge;
	print "a:",a;
	print "b:",b;
	print "Points:",Points;
	*/
	if Edge[1] lt Edge[2] then
		PP := Remove(Remove(Points,Edge[1]),Edge[2]-1);
	else
		PP := Remove(Remove(Points,Edge[2]),Edge[1]-1);
	end if;
	// Make vector of centers of circumcircles
	PPP := [ (2*PP[k]-b-a)/(b-a) : k in [1..#PP] ];
	
	// Factors due to change of variable
	Fact1 := (b-a)/2;
	Fact2 := (b+a)/(b-a);
	/*
	print "PP:",PP;
	print "PPP:",PPP;
	print "Fact1:",Fact1;
	print "Fact2:",Fact2;
	print "DFF:",DFF;
	*/
	vprint SE,2 : "########### Next Integral ###################";
	Integrals0 := [ CC_0 : j in [1..g] ];
	for t in [1..#Abscissas] do
		x := Abscissas[t];

		// Analytic continuation
		z := RS_NRootAC(x,PPP,Zeta,N);
		//print "z:",z;
		//print "Arg(z):",Arg(z);
		Enum := (x + Fact2);
		for j in [1..g] do
			w := DFF[j];
			//Factor := Fact1^((N*(w[1]+1)-d*w[2])/N);
			//print "((w[1]+1)-(d*w[2]/N)):",((w[1]+1)-(d*w[2]/N));
			//print "z:",z;
			Denom := z^w[2];
			//print "(2*(N-w[2])/N):",(2*(N-w[2])/N);
			dx := Weights[t][1] * Weights[t][2]^(2-(2*w[2]/N));
			Integrals0[j] +:= (Enum^w[1]/Denom) * dx;
		end for;
	end for;
	Integrals := []; // g x N array of integrals
	for j in [1..g] do
		w := DFF[j];
		Factor := Fact1^((w[1]+1)-(d*w[2]/N));
		Integrals0[j] *:= (Exp(-Pi*i*w[2]/N)-Exp(Pi*i*w[2]/N)) * StepLength * Factor;
		Append(~Integrals,[ Zeta^(l*w[2]) * Integrals0[j] : l in [0..N-2]]);
	end for;
		
	return Integrals;
end function;


intrinsic RS_SEPM( f::RngMPolElt : Prec := -1, Small := true ) -> Matrix
{ Computes period matrices associated to the superelliptic curve defined by f to precision Prec }
	// Precision management
	GlobalPrec := RS_GetGlobalPrecision();
	if Prec lt 20 then
		Prec := GlobalPrec;
	end if;
	vprint SE,1 : "Prescribed precision:",Prec;
	

	// Branch points
	Points := Sort(RS_DiscriminantPoints(f,Parent(f).2:Prec:=Prec),RS_CompareFldComElt);
	vprint SE,3 : "Discriminant points:",Points;
	vprint SEP,1 : "Precision of discriminant points:",Precision(Points[1]);;

	// Degrees
	N := Degree(f,2); d := Degree(f,1);
	//assert &and[N ge 2, d ge 2, N*d ge 6]; // No Genus 0 curves
	assert &and[N ge 2, d ge 3];

	// Complex fields and constants
	C_<i> := ComplexField(Prec); C_20<i> := RS_GetGlobalComplexField_20(); RS_SetGlobalZeta(N); Zeta := RS_GetGlobalZeta(); 
	CC<i> := RS_GetGlobalComplexField_Max(); CC_0 := Zero(CC); Pi := RS_GetGlobalPi(); 
	vprint SE,1 : "Zeta:",C_20!Zeta;
	vprint SEP,1 : "Precision Zeta:",Precision(Zeta);
	vprint SEP,1 : "Precision CC_0:",Precision(CC_0);
	vprint SEP,1 : "Precision Pi:",Precision(Pi);


	// Low precision branch points
	LowPrecPoints := [ C_20!Pt : Pt in Points ];

	// Maximal spanning tree w.r.t holomorphicy
	vprint SE,1 : "Constructing maximal spanning tree...";
	t := Cputime();
	MST_Edges, MST_Tau := RS_SEMST(LowPrecPoints);
	Cputime(t);
	vprint SE,3 : "MST_Edges:",MST_Edges;
	vprint SE,2 : "#MST_Edges:",#MST_Edges;
	vprint SE,2 : "MST_Tau:",MST_Tau;

	"Prec!!!!:",RS_GetGlobalPrecision();

	// Holomorphic differentials
	vprint SE,1 : "Computing holomorphic differentials...";
	t := Cputime();
	DFF := [];
	for i in [0..d-2] do
		for j in [1..N-1] do
			//if (-N*(i+1) + j*d - Gcd(N,d) ge 0) and (N*(i+1) - j - 1 ge 0) then
			if (-N*(i+1) + j*d - Gcd(N,d) ge 0) then
				Append(~DFF,<i,j>);
			end if;
		end for;
	end for;
	vprint SE,2 : "Holomorphic differentials:",DFF;
	Cputime(t);

	// Genus
	g := #DFF;
	assert 2*g eq (N-1)*(d-1) - Gcd(N,d) + 1; // Check Riemann-Hurwitz
	vprint SE,1 : "Genus:",g;

	// Integration parameters
	vprint SE,1 : "Computing Integration parameters...";
	t := Cputime();
	Abscissas, Weights, StepLength := RS_SEIntegrationParameters(Points,MST_Edges,MST_Tau,N);
	vprint SE,2 : "#Abscissas:",#Abscissas;
	vprint SE,3 : "Abscissas:",Abscissas;
	vprint SE,2 : "#Weights:",#Weights;
	vprint SE,3 : "Weights:",Weights;
	vprint SE,2 : "StepLength:",C_20!StepLength;
	vprint SEP,1 : "Precision Abscissas:",Precision(Abscissas[1]);
	vprint SEP,1 : "Precision Weights:",Precision(Weights[1][1]);
	vprint SEP,1 : "Precision StepLength:",Precision(StepLength);
	Cputime(t);
"Prec!!!!:",RS_GetGlobalPrecision();
	// Integrals
	vprint SE,1 : "Integrating...";
	t := Cputime();
	Integrals := [];
	for k in [1..d-1] do
		vprint SE,2 : "Integrating edge nr.",k;
		I := RS_SEIntegrate( MST_Edges[k],Points,DFF,Abscissas,Weights,StepLength,d,N );
		Append(~Integrals,I);
	end for;
	vprint SEP,1 : "Precision Integrals:",Precision(Integrals[Random([1..d-1])][Random([1..g])]);
	Cputime(t);

	// Matrix of periods
	//PM_ := ZeroMatrix(CC,g,2*g);
	//PM_ := ZeroMatrix(CC,g,(d-1)*(N-1));
	PM_ := ZeroMatrix(CC,(d-1)*(N-1),g);
	for k in [1..d-1] do
		for l in [1..N-1] do
		//for l in [1..N] do
			Ind2 := (N-1)*(k-1) + l;
			//Ind2 := N*(k-1) + l;
			for j in [1..g] do
				//PM_[j][Ind2] := Integrals[k][j][l];
				PM_[Ind2][j] := Integrals[k][j][l];
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
			skd := RS_SEIntersection(MST_Edges[j],MST_Edges[l],Points,N);
			vprint SE,2: "<s,k,d>-triple of edges",MST_Edges[j],"and",MST_Edges[l],":",skd;
			skd_Matrix[j][l] := <skd[1],skd[2] mod N,skd[3]>;
			skd_Matrix[l][j] := <-skd[1],skd[2] mod N,skd[3]>;
		end for;
	end for;
	Cputime(t);
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
	RealIntMatrix := BlockMatrix(d-1,d-1,Blocks);
	
	vprint SE,3: "(Real) Intersection matrix:",RealIntMatrix;
	vprint SE,2: "Rank(IntersectionMatrix):",Rank(RealIntMatrix);
	assert Rank(RealIntMatrix) eq 2*g;	
	
	// Symplectic reduction of intersection matrix
	vprint SE,1 : "Performing symplectic reduction...";
	t := Cputime();
	//CF, Alpha := RS_SymplecticBasis(RealIntMatrix);
	CF, Alpha := RS_SymplecticBasis(Transpose(RealIntMatrix));
	Cputime(t);
	vprint SE,2: "Alpha:",Alpha;
	vprint SE,3: "CF:",CF;
	
	//Alpha_CC := ChangeRing(Transpose(Alpha),CC);
	Alpha_CC := ChangeRing(Alpha,CC);
	vprint SE,1 : "Matrix multiplication 1...";
	t := Cputime();
	//PMAPMB := PM_ * Alpha_CC;
	PMAPMB := Alpha_CC * PM_;
	Cputime(t);
	
	//PM_A := ColumnSubmatrixRange(PMAPMB,1,g);
	//PM_B := ColumnSubmatrixRange(PMAPMB,g+1,2*g);
	PM_A := RowSubmatrixRange(PMAPMB,1,g);
	PM_B := RowSubmatrixRange(PMAPMB,g+1,2*g);
	

	//vprint SE,3 : "Dependent columns:",ColumnSubmatrixRange(PMAPMB,2*g+1,Nrows(Alpha));
	vprint SE,3 : "Dependent rows:",RowSubmatrixRange(PMAPMB,2*g+1,Nrows(Alpha));
	vprint SE,3 : "PM_A:",PM_A;
	vprint SE,3 : "PM_B:",PM_B;
	//vprint SE,2 : "Det(PM_A):",Determinant(PM_A);
	//vprint SE,2 : "Det(PM_B):",Determinant(PM_B);
	//vprint SE,3 : "PM_B^-1:",PM_B^-1;
	

	if Small eq false then
		// Compute big period matrix
		PM := HorizontalJoin(PM_B,PM_A) + ZeroMatrix(C_,g,2*g); // Dimensions???
		vprint SE,1 : "Period matrix:";
		RS_SetGlobalPrecision(GlobalPrec);
		return PM;
	end if;


	vprint SE,1 : "Matrix inversion...";
	t := Cputime();
	//PM_BInv := PM_B^(-1);
	PM_AInv := PM_A^(-1);
	Cputime(t);
	vprint SE,1 : "Matrix multiplication 2...";
	t := Cputime();
	//PM := PM_BInv * PM_A;
	//PM := PM_A * PM_BInv;
	PM := PM_B * PM_AInv;
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
		print "f:",f;
		error Error("Period matrix not symmetric: Computation failed.");
	end if;	
		
	// Testing positive definiteness of the imaginary part of the period matrix
	vprint SE,1 : "Testing positive definite...";
	t := Cputime();
	PM_Im := ZeroMatrix(RealField(Prec),g,g);
	for j in [1..g] do
		for k in [j..g] do
			PM_Im[j][k] := Real(Im(PM[j][k]));
			PM_Im[k][j] := PM_Im[j][k];
		end for;
	end for;
	assert IsPositiveDefinite(PM_Im);
	Cputime(t);

	RS_SetGlobalPrecision(GlobalPrec);

	return PM;

end intrinsic;



intrinsic RS_SuperellipticDifferentials( f::RngMPolElt ) -> SeqEnum[DiffFunElt]
{ Returns a basis of holomorphic differentials for superelliptic curves }
	N := Degree(f,2);
	d := Degree(f,1);
	DFF := [];
	for i in [0..d-2] do
		for j in [1..N-1] do
			if (-N*(i+1) + j*d - Gcd(N,d) ge 0) then
				Append(~DFF,<i,j>);
			end if;
		end for;
	end for;
	return DFF;
end intrinsic;

function RS_SEIntegrationPoints( h, NrPoints, Lambda )
// Computes integration parameters for double exponential integration on the interval [-1,1]
	// Precision of integration parameters only depens on Precision(h)
	eh := Exp(h); // e^h
	eh_inv := 1/eh; // e^-h
	ekh := 1; // e^0h
	ekh_inv := 1; // e^-0h
	Abscissas_ := [];
	Weights_ := [];
	for k in [1..NrPoints] do
		ekh *:= eh; // e^kh
		ekh_inv *:= eh_inv; // e^-kh
     		sh := (ekh-ekh_inv)/2; // sinh(kh) = (1/2) * (e^kh - e^-kh)
		ch := (ekh+ekh_inv)/2 ; // cosh(kh) = (e^kh + e^-kh)
      		//ch2 := (ekh+ekh_inv); //  2*cosh(kh) = (e^kh + e^-kh)
      		esh := Exp(Lambda*sh); // e^(Lambda*sinh(kh))
      		esh_inv := 1/esh; // e^-(Lambda*sinh(kh))
		chsh_inv := 2/(esh+esh_inv); // 1/cosh(Lambda*sinh(kh)))
      		//chsh2_inv := 1/(esh+esh_inv); // 1/(2*cosh(Lambda*sinh(kh)))
		shsh := (esh-esh_inv)/2; // sinh(Lambda*sinh(kh))
      		//shsh2 := esh-esh_inv; // 2*sinh(Lambda*sinh(kh))
		thsh := shsh*chsh_inv; // tanh(Lambda*sinh(kh)) =  sinh(Lambda*sinh(kh)) / cosh(Lambda*sinh(kh))
      		//thsh := shsh2*chsh2_inv; // tanh(Lambda*sinh(kh)) =  2*sinh(Lambda*sinh(kh)) / 2*cosh(Lambda*sinh(kh))
		Append(~Abscissas_,thsh);
		Append(~Weights_,<ch,chsh_inv>); // = <cosh(kh),1/cosh(Lambda*sinh(kh))>
		//Append(~Weights_,2*ch2*chsh2_inv^2); // = 2*2*cosh(kh) / 4*(cosh(Lambda*sinh(kh)))^2
     	end for;

	Abscissas := Reverse([-zk : zk in Abscissas_]) cat [0] cat Abscissas_;
	Weights := Reverse(Weights_) cat [<1,1>] cat Weights_;

  	return Abscissas, Weights, h*Lambda;
	//return Abscissas, Weights, h;
end function;

function RS_Dist_1( P )
// Strange distance function?!?
	xP := Abs(Re(P));
	i := Parent(P).1;
	if xP gt 1 then
		return Abs(xP-1+i*Im(P));
	else
		return Abs(Im(P));
	end if;
end function;

function RS_M1Bound(PPP,N)
// Compute the bound M1 on the integrand
	B1 := 1/(&*[ RS_Dist_1(P) : P in PPP ])^(1/N);
	//print "B1:",B1;
	if B1 gt 1 then
		return B1^(N-1);
	else
		return B1;
	end if;
end function;


function RS_tOpt(D,M1)
// Some value
	Pi := RS_GetGlobalPi();
	Lambda := Pi/2;
	return Argsinh((D+Log(4*M1))/Lambda);
end function;


function RS_Dist_thsh( P, Tau )
// Strange distance function?!?
	Pi := RS_GetGlobalPi();
	Lambda := Pi/2;
	i := Parent(P).1;
	P := Abs(Re(P)) + i*Abs(Im(P));
	x0 := 0; x1 := Argcosh(Pi/(2*Lambda*Sin(Tau))); // s.t. \Lambda Cosh(x)Sin(\Tau)= Pi/2;
	Phi := function(x)
		return Tanh(Lambda*Sinh(x+i*Tau));
	end function;
	ArgPhiP := function(x)
		return Arg(Cosh(x+i*Tau))-2*Arg(Cosh(Lambda*Sinh(x+i*Tau)));
	end function;
	ArgPhi := function(x)
		return Arg(P - Phi(x));
	end function;
	// x := Solve(x := x0, x1, ArgPhi(x)-ArgPhiP(x)-Pi/2);
	x := x0;
	Val := ArgPhi(x)-ArgPhiP(x)-Pi/2;
	//print "x0:",x;
	//print "x1:",x1;
	//print "Val:",Val;
	n := 10;
	for t in [0..n] do
		x_new := (t/n)*x1 + (1-(t/n))*x0;
		Val_new := ArgPhi(x_new)-ArgPhiP(x_new)-Pi/2;
		if Abs(Val_new) lt Abs(Val) then
			x := x_new;
			Val := Val_new;
			
		end if;
	end for;
	//print "x_new:",x_new;
	//print "Val_new:",Val_new;
	return Abs(P-Phi(x));
end function;


function RS_M2Bound(PPP,Tau,N)
// Compute the bound M2 on the integrand
	B2 := 1/(&*[ RS_Dist_thsh(P,Tau) : P in PPP ])^(1/N);
	//print "B2:",B2;
	if B2 gt 1 then
		return B2^(N-1);
	else
		return B2;
	end if;
end function;


function RS_PhiBound(Tau,N)
// .. 
	Pi := RS_GetGlobalPi();
	Lambda := Pi/2;
	Alpha := (N-1)/N;
	Y_Tau := Sqrt(Pi/2*Lambda*Sin(Tau));
	X_Tau := Sqrt(Y_Tau^2 - Lambda^2 * Sin(Tau)^2 / Tan(Tau));
	return ((1+2*Alpha)/(2*Alpha*Cos(Y_Tau)^(2*Alpha))+1/(2*Alpha*Sinh(X_Tau)^(2*Alpha)));
	// something different...
end function;


function RS_hOpt(D,Tau,M2,N)
	Pi := RS_GetGlobalPi();
	Lambda := Pi/2;
	return 2*Pi*Tau/Log(1+2*M2*RS_PhiBound(Tau,N)*Exp(D));
end function;


intrinsic RS_SEIntegrationParameters( Points::SeqEnum[FldComElt], MST::SeqEnum[Tup], Tau::FldReElt, N::RngIntElt ) -> SeqEnum[FldReElt], SeqEnum[FldReElt], FldReElt
{ Computes integration parameters for double exponential integration on the interval [-1,1] }
	assert N ge 2;

	// Get complex field for printing
	C_20<i> := RS_GetGlobalComplexField_20();
	// Get complex field of maximal precision
	C<i> := RS_GetGlobalComplexField_Max(); 
	// Get Pi
	Pi := RS_GetGlobalPi();
	// Compute D
	D := Real(C!(Precision(RS_GetGlobalComplexField_Comp()) * Log(10)));"PD:",Precision(D);
	"Prec!222!!:",RS_GetGlobalPrecision_Comp();
	// Check Tau
	if Tau ge Pi/2 then
		//Tau *:= (19/20);
		error Error("OKAY?!?!?");
	end if;
	assert &and[Tau lt Pi/2, Tau gt 0];

	// Coerce into larger field make it precise
	Tau *:= (19/20);
	Tau := C!Tau;
	vprint Int,1 : "Tau:",C_20!Tau;

	/*
	// Compute bounds M1,M2
	M1 := 0; M2 := 0;
	for Edge in MST do
		if Edge[1] lt Edge[2] then
			PP := Remove(Remove(Points,Edge[1]),Edge[2]-1);
		else
			PP := Remove(Remove(Points,Edge[2]),Edge[1]-1);
		end if;
		// Make vector of centers of circumcircles
		PPP := [ (2*PP[k]-Points[Edge[2]]-Points[Edge[1]])/(Points[Edge[2]]-Points[Edge[1]]) : k in [1..#PP] ];
		M1 := Max(M1,RS_M1Bound(PPP,N));
		M2 := Max(M2,RS_M2Bound(PPP,Tau,N));
	end for;
	M1 := Ceiling(M1); M2 := Ceiling(M2);
	vprint Int,1 : "M1:",M1;
	vprint Int,2 : "M2:",M2;

	h := RS_hOpt(D,Tau,M2,N);
	n := Ceiling(RS_tOpt(D,M1)/h);
	*/
	
	// Strange parameter lambda -> Explained in Pascal's thesis
	Lambda := Pi/2;

	// Other parameters
	
	I := 100000; // What is I?
	//Alpha := (N-1)/N;
	Alpha := 1/N;
	f_Max := 1000000;
	Y_Tau := (Pi/2) * Sin(Tau)^(7/8);
	X_Tau := (1/Tan(Tau)) * Sqrt( Y_Tau^2 - ( Pi*Sin(Tau) / 2)^2 );
	M_Tau := (2 * f_Max / Cos(Tau) ) * ( Tanh(X_Tau)/Cos(Y_Tau)^2  +  1/Tanh(X_Tau));

	// Compute step length h
	h := Real(C!(2*Pi * Tau ) / ( D + Log( ((16*X_Tau * M_Tau * I)/Cos(Tau)) + 1) ));
	vprint Int,2 : "f_Max:",C_20!f_Max;
	vprint Int,2 : "X_Tau:",C_20!X_Tau;
	vprint Int,2 : "Y_Max:",C_20!Y_Tau;
	vprint Int,2 : "M_Max:",C_20!M_Tau;
	// Compute #integration points = 2*n+1
	print "...:",Argsinh( D+ Log((4^(1-Alpha) * M_Tau)/(2-2*Alpha))/(Lambda*(2-2*Alpha)) ) / h;
	n := Ceiling(Argsinh( D+ Log((4^(1-Alpha) * M_Tau)/(2-2*Alpha))/(Lambda*(2-2*Alpha)) ) / h );
	
	"PH:",Precision(h);
	"PH:",Precision(Pi);
	"PH:",Precision(Tau);
	"PH:",Precision(D);
	vprint Int,1 : "h:",C!h;
	vprint Int,1 : "n:",n;

	// Compute integration points
	Abscissas, Weights, StepLength := RS_SEIntegrationPoints(h,n,Lambda);
	
	return Abscissas, Weights, StepLength;

end intrinsic;


// ################################################################  !!! FROM HERE !!! ################################################################
// ################################################################  !!! OLD  CODE !!! ################################################################

function RS_SENewtonIteration(f,df,z,k )
	if k gt 0 then
		return RS_SENewtonIteration(f,df,z - (Evaluate(f,z)/Evaluate(df,z)) ,k-1);
	end if;
	return z;
end function;


intrinsic RS_SENewtonMethod( f::RngUPolElt, z::FldComElt ) -> FldComElt
{ Approximate a zero of f=df[1] using newton's method starting with ApproxValues }
	
	// Compute upper bound on the size of the roots of f
	Bound := RS_UpperBoundRoots(f);
	vprint RA,1 :"Bound on roots:",Bound;

	// Compute required # of extra digits
	ExtraDigits := RS_ExtraDigits(Bound);
	
	// Adjusted complex field and polynomial ring for calculations
	C_<i> := ComplexField(Precision(RS_GetGlobalComplexField_Comp())+ExtraDigits+10); C_y<y> := PolynomialRing(C_); R_1 := Real(One(C_));

	// Coercion of root approximation data
	if Precision(z) lt Precision(C_) then
		z := C_!z;
	end if;
	
	if Precision(BaseRing(f)) lt Precision(C_) then
		f := C_y!f;
	end if;
	
	// Start root approximation
	BetaValues := [];
	N := Degree(f);
	df := [f];
	for k in [1..N] do
		Append(~df,Derivative(df[#df]));
	end for;

	df := [f];
	for k in [1..Degree(f)] do
		Append(~df,Derivative(df[#df]));
	end for;

	// Check if z is an approximate zero of f(x_t2,y)
	Okay, Beta_z_j := RS_AlphaTest(df,z);

	if not Okay then
		return 0;
	end if;
	
	return RS_SENewtonIteration( df[1],df[2],z,RS_GetRAIterations() );
end intrinsic




function RS_SEAC(f,Gamma,t1,t2,y_t1)
// Analytically continue the fiber above Gamma(t1) to the fiber above Gamma(t2)
	// Evaluate the path Gamma at t1 and t2
	x_t2 := RS_PathEvaluate(Gamma,(1/2)*(t2+1));
	
	// Evaluate y-derivatives of f at x_t2
	f_x_t2 := RS_MPolEval(f,x_t2);

	// Divide by leading coefficient
	f_x_t2 *:= 1/LeadingCoefficient(f_x_t2);

	
	vprint AC,3 : "t1:",t1;
	//vprint AC,2: "y_t1:",y_t1;
	vprint AC,3 : "t2:",t2;
	
	y_t2 := RS_SENewtonMethod(f_x_t2,y_t1);

	if y_t2 ne 0 then
		return y_t2;
	else
		//vprint AC,2: "Deeper recursion: Not approximate or no convergence to different solutions.";
		y_t2_1 := RS_SEAC(f,Gamma,t1,(t1+t2)/2,y_t1);
		y_t2_2 := RS_SEAC(f,Gamma,(t1+t2)/2,t2,y_t2_1);
		return y_t2_2;
	end if;
end function;

function RS_SEACRecursion(f,Gamma,t1,t2,Fiber_t1)
// Analytically continue the fiber above Gamma(t1) to the fiber above Gamma(t2)

	vprint PAC,1 : "Precision(t1):",Precision(t1);
	vprint PAC,1 : "Precision(t2):",Precision(t2);

	// Evaluate the path Gamma at t1 and t2
	//x_t1 := RS_PathEvaluate(Gamma,t1);
	//x_t2 := RS_PathEvaluate(Gamma,t2);
	x_t2 := RS_PathEvaluate(Gamma,(1/2)*(t2+1));
	vprint PAC,1 : "Precision(x_t2):",Precision(x_t2);

	// Evaluate f(x,y) at x = x_t2
	//f_x_t2 := UnivariatePolynomial(Evaluate(f,1,x_t2));
	//f_x_t2 := RS_MPolEval(f,x_t2);
	f_x_t2 := Evaluate(f,x_t2);
	vprint PAC,1 : "Precision(f_x_t2):",Precision(BaseRing(Parent(f_x_t2)));

	vprint AC,2 : "Gamma:",Gamma;
	vprint AC,2 : "x_t2:",x_t2;
	vprint AC,3 :"f_x_t2:",f_x_t2;

	// Divide by leading coefficient
	f_x_t2 *:= 1/LeadingCoefficient(f_x_t2);
	/*if LeadingCoefficient(f_x_t2) ne 1 then
		RS_SetTmp(<f_x_t2,Fiber_t1,t1,t2>);
		return "bla1";
	end if;*/
	vprint PAC,1 : "Precision fxt2:",Precision(LeadingCoefficient(f_x_t2));	

	vprint AC,3 : "t1:",t1;
	vprint AC,3 : "t2:",t2;
	vprint AC,3 : "Difference:", t2-t1;
	if t2-t1 eq 0 then
		RS_SetTmp(<f_x_t2,Fiber_t1,t1,t2>);
		//return "bla";
	end if;

	S := RS_GetRAMethod();
	vprint PAC,1 : "Precision(RS_GetGlobalComplexField()):",Precision(RS_GetGlobalComplexField());
	vprint PAC,1 : "Precision(RS_GetGlobalComplexField_Comp()):",Precision(RS_GetGlobalComplexField_Comp());
	vprint PAC,1 : "Precision(RS_GetGlobalComplexField_Max()):",Precision(RS_GetGlobalComplexField_Max());
	
	if S eq "DKPEB" then
		Fiber_t2 := RS_DurandKernerPEB(f_x_t2,Fiber_t1,Precision(RS_GetGlobalComplexField_Comp()));
	elif S eq "BSWPEB" then
		Fiber_t2 := RS_BoerschSupanWeierstrassPEB(f_x_t2,Fiber_t1,Precision(RS_GetGlobalComplexField_Comp()));
	elif S eq "DK" then
		Fiber_t2 := RS_DurandKerner(f_x_t2,Fiber_t1);
	elif S eq "BSW" then
		Fiber_t2 := RS_BoerschSupanWeierstrass(f_x_t2,Fiber_t1);
	elif S eq "EAN" then
		Fiber_t2 := RS_EhrlichAberthNewton(f_x_t2,Fiber_t1);
	elif S eq "Newton" then
		Fiber_t2 := RS_NewtonMethod(f_x_t2,Fiber_t1);
	end if;

	if #Fiber_t2 gt 0 then
		return Fiber_t2;
	else
		vprint AC,2: "Deeper recursion: Not approximate or no convergence to different solutions.";
		Fiber_t2_1 := RS_ACRecursion(f,Gamma,t1,(t1+t2)/2,Fiber_t1);
		vprint AC,2: "Fiber_t2_1:",Fiber_t2_1;
		Fiber_t2_2 := RS_ACRecursion(f,Gamma,(t1+t2)/2,t2,Fiber_t2_1);
		vprint AC,2: "Fiber_t2_2:",Fiber_t2_2;
		return Fiber_t2_2;
	end if;
end function;

intrinsic RS_SEPM2( f::RngMPolElt : Prec := -1, Basepoint := "Opt", Big := false, PathMethod := RS_GetPathMethod(), RAMethod := "Newton", IntMethod := "DE" ) -> Mtrx
{ Computes the period matrix of the Riemann surface associated to the superelliptic curve defined by the polynomial f(x,y) = y^N - g(x) }

	vprint PM,1 : "Defining polynomial: ",f;

	// Checking for special cases
	f_y := UnivariatePolynomial(Evaluate(f,1,1));
	f_y *:= 1/LeadingCoefficient(f_y);
	if f_y - LeadingTerm(f_y) in BaseRing(f_y) then
		f_x := UnivariatePolynomial(Evaluate(f,2,0)); 
		if Gcd(f_x,Derivative(f_x)) eq One(Parent(f_x)) then
			if Degree(f_y) eq 2 then
				if Degree(f_x) eq 3 then
					print "Special case: Curve is elliptic!";
				else
					print "Special case: Curve is hyperelliptic!";
				end if;
			else
				print "Special case: Curve is superelliptic!";
			end if;
		else
			print "Curve is not superelliptic!";
			return 0;
		end if;
	else
		print "Curve is not superelliptic!";
		return 0;
	end if;

	// Computing the a basis of homolomorphic differentials
	
	FF<x,y> := FunctionField(f);
	//DFF := BasisOfDifferentialsFirstKind(FF);
	DFF := RS_SuperellipticDifferentials(f);
	vprint PM,3 : "Basis of homolormphic 1-forms: ",DFF;

	// Genus
	g := #DFF;
	vprint PM,1 : "Genus of the Riemann surface:",g;


	// Dealing with genus 0 curves
	if g eq 0 then
		print "Genus of the Riemann surface is zero!";
		return [];
	end if;	


	// Representatives of holomorphic differentials
	//vprint PM,1 : "Computing representatives:";
	//time DFF_Reps := [ RationalFunction( (DFF[j]/ Differential(FF!x)),Rationals() ) : j in [1..g] ];
	//DFF_Reps := [ RationalFunction(DFF[j],Rationals() ) : j in [1..g] ];
	DFF_Reps := DFF;
	vprint PM,2 : "DFF_Reps:",DFF_Reps;	

	// Discriminant Points
	DiscriminantPoints := RS_DiscriminantPoints(f,Parent(f).2);

	// #Sheets
	N := Degree(f,2);
	Sym := Sym(N);
	Id := Id(Sym);
	vprint PM,1 : "#Sheets:",N;


	 N_x := Degree(f,2);
	
	// Precision management
	GlobalPrec := RS_GetGlobalPrecision();
	if Prec lt 20 then
		Prec := GlobalPrec;
	end if;

	// Adjust Prec "somehow" (yet to be done!)
	// Compute maximal loss of Prec after integrating
	NumberOfCycles := 2*g + N -1;
	Adjusted := 10;

	// Set precision
	RS_SetGlobalPrecision( Prec : Adjusted := Adjusted );

	// Complex field for return
	C_<i> := RS_GetGlobalComplexField();
	vprint PM,1 : "Prescribed precision:",Precision(C_);

	// Complex field used for computations
	C<i> := RS_GetGlobalComplexField_Comp(); Cy<y> := PolynomialRing(C); Cxy<x,y> := PolynomialRing(C,2);
	vprint PM,1 : "Adjusted precision:",Precision(C); 

	// Complex field for printing
	C_20<i> := RS_GetGlobalComplexField_20();

	// Complex field of maximal precision
	CC<i> := RS_GetGlobalComplexField_Max(); CC_0 := Zero(CC); CC_1 := One(CC); RR_0 := Real(CC_0); RR_1 := Real(CC_1); MaxPrec := Precision(CC); 
	CCy<y> := PolynomialRing(CC); CCxy<x,y> := PolynomialRing(CC,2); f := CCxy!f;
	vprint PM,1 : "Maximal precision:",MaxPrec;

	// Error control
	Err := RS_GetGlobalError();
	vprint PM,1 : "GlobalError:",C_20!Err;

	// Compute roots of unity
	Zeta := RS_ReducePrecision(Exp(2*RS_GetGlobalPi()*CC.1/N),5); ZetaVec := [ Zeta^k : k in [0..N-1] ];
	vprint PM,1 : "Primitve root of unity:",Zeta;
	vprint PM,3 : "generates:",ZetaVec;
	
	// Compute fundamental group w.r.t PathMethod
	vprint PM,1 :"Computing fundamental group using path method:",PathMethod;
	if IntMethod eq "GQ" then
		RS_SetMaxSafeRadius(1/1);
	elif IntMethod eq "DE" then
		RS_SetMaxSafeRadius(1/500);
	end if;
	Basepoint, OrdDiscPts, PathPieces, IndexPathLists := RS_PathPieces( DiscriminantPoints : Basepoint := Basepoint, PathMethod := PathMethod );
	
	assert Precision(Zeta) eq MaxPrec - 5;
	assert Precision(PathPieces[1]`EndPt) eq MaxPrec - 10;
	assert Precision(Basepoint) eq MaxPrec - 5;
	assert Precision(DiscriminantPoints[1]) eq MaxPrec - 5;

	// Compute length of paths
	PathLengths := [ RS_PathLength(PathPieces[j]) : j in [1..#PathPieces] ];
	TotalPathLength := &+PathLengths;
	MinPathLength := Min(PathLengths);
	MaxPathLength := Max(PathLengths);

	vprint PM,1 : "Total pathlength:",C_20!TotalPathLength;	
	vprint PM,1 : "Minimal path length:",C_20!MinPathLength;
	vprint PM,1 : "Maximal path length:",C_20!MaxPathLength;

	vprint PM,1 : "Basepoint:",C_20!Basepoint;
	vprint PM,1 : "#Discriminant points:",#OrdDiscPts;
	vprint PM,3 : "Discriminant points:",OrdDiscPts;
	vprint PM,1 : "#Integrals:",#PathPieces;
	vprint PM,3 : "PathPieces:",PathPieces;

	// Compute Integration parameters
	assert IntMethod in ["DE","GQ"];
	vprint PM,1 :"Using integration method:",IntMethod;
	if IntMethod eq "DE" then
		// Compute Tau for integration
		LowPrecOrdDiscPts := [ C_20!P : P in OrdDiscPts ];
		RS_TauPaths( LowPrecOrdDiscPts , PathPieces );
		vprint PM,3 : [ <Gamma`Tau,Gamma`Type> : Gamma in PathPieces ];
		MinTau, Ind := Min( [ Gamma`Tau : Gamma in PathPieces ] );
		vprint PM,1 : "Minimal Tau realized for path:",PathPieces[Ind],"of length:",C_20!RS_PathLength(PathPieces[Ind]);
		RS_SetIntTau((19/20)*Real(CC!MinTau));
		//RS_SetIntTau( RS_GetGlobalPi()/10);
		"#######################################################";
		vprint PM,1 : "Tau(integration):",C_20!RS_GetIntTau();
		"#######################################################";
		assert Precision(RS_GetIntTau()) eq Precision(CC_0);
		vprint PM,1 :"Computing integration parameters:";
		time Abscissas, Weights, StepLength := RS_IntegrationParameters(RS_GetIntTau());
	elif IntMethod eq "GQ" then
		vprint PM,1 :"Computing integration parameters:";
		Steps := Ceiling((7/2)*Prec*Log(Prec));
		//Steps := Ceiling((7/2)*Prec); //Steps := 1000;
		StepLength := 1/Steps;
		time Abscissas, Weights := RS_GaussQuadratureParameters(Steps); 
	end if;
	
	vprint PM,1 : "Integrating in ",#Abscissas,"Steps and Steplength",C_20!StepLength;
	vprint PM,3 : "Abscissas:",Abscissas; 	
	vprint PM,3 : "Weights:",Weights;
	assert Precision(Abscissas[1]) le MaxPrec;
	assert Precision(Weights[1]) eq Precision(Abscissas[1]);
	assert Precision(Abscissas[1]) ge Precision(C);


	// Compute optimal constant and number of iterations for root approximation
	/*
	if N eq 2 then
		RS_SetRAMethod("DKPEB",N);
	else
		RS_SetRAMethod(RAMethod,N);
	end if;
	*/
	RS_SetRAMethod(RAMethod,N);
	vprint PM,1 : "Using root approximation method:",RS_GetRAMethod();

	// Analytically continue each path along each differential form
	for l in [1..#PathPieces] do
		
		Gamma := PathPieces[l];
		
		vprint PM,2 : "#####################################################################";
		vprint PM,2 : "Integrating path (#",l,"):",Gamma;

		y_0 := RS_Fiber(f,Gamma`StartPt);
		/*
		Gamma_y_Values := [];
		
		//y_j := RS_ACRecursion(f,Gamma,R_0,Abscissas[1],y_0); 
		y_j := RS_ACRecursion(f,Gamma,-R_1,Abscissas[1],y_0);

		Append(~Gamma_y_Values,y_j);

		//print "Prec(Absc):",Precision(Abscissas[1]);
		for j in [2..#Abscissas] do
			y_j := RS_ACRecursion(f,Gamma,Abscissas[j-1],Abscissas[j],y_j);
			Append(~Gamma_y_Values,y_j);
		end for;
		

		//print "Prec y_j:",Precision(y_j);	

		y_j := RS_ACRecursion(f,Gamma,Abscissas[#Abscissas],R_1,y_j);
		*/
		
		Gamma_y_Values := [];

		y_0 := RS_Fiber(f,Gamma`StartPt);
		y_0_RoU := [ y_0[1]*ZetaVec[j] : j in [1..N] ];

		_, Perm_RoU := RS_Permutation(y_0,y_0_RoU,RS_GetGlobalError()); // Make this error = 0
		
		ZetaVec := [ ZetaVec[j^Inverse(Perm_RoU) ] : j in [1..N] ];
		//y_0_RoU := [ y_0[1]*ZetaVec[j] : j in [1..N] ];
		//_, Perm_RoU := RS_Permutation(y_0,y_0_RoU,RS_GetGlobalError()); // Make this error = 0
		//assert Perm_RoU eq Id;

		// Continue to first absciss
		y_j_ := RS_SEAC(f,Gamma,-RR_1,Abscissas[1],y_0[1]); 
		y_j := [ y_j_*ZetaVec[j] : j in [1..N] ];

		Append(~Gamma_y_Values,y_j);

		for j in [2..#Abscissas] do
			y_j_ := RS_SEAC(f,Gamma,Abscissas[j-1],Abscissas[j],y_j[1]);
			y_j := [ y_j_*ZetaVec[j] : j in [1..N] ];
			Append(~Gamma_y_Values,y_j);
		end for;	
		
		// Continue to next start point
		y_j_ := RS_SEAC(f,Gamma,Abscissas[#Abscissas],RR_1,y_j[1]);
		y_j := [ y_j_*ZetaVec[j] : j in [1..N] ];
		


		// Compute the local monodromy of the path...
		Ok, Sigma := RS_Permutation(RS_Fiber(f,Gamma`EndPt),y_j,RS_GetGlobalError());

		if not Ok then
			print "DefiningPolynomial:",f;
			RS_SetGlobalPrecision(GlobalPrec);
			error Error("Sheets do not match after analytic continuation.");
		end if;	
		
		vprint PM,2 : "Sigma:",Sigma;

		if assigned Gamma`Permutation then
			if Gamma`Permutation ne Sigma then
				print "DefiningPolynomial:",f;
				RS_SetGlobalPrecision(GlobalPrec);
				error Error("Different local monodromies computed while integrating.");
			end if;
		else
			// Assign the permutation
			Gamma`Permutation := Sigma;
		end if;	

		// Now integrate
		PathDiffIntegrals := [];
		for Rep in DFF_Reps do
			PathDiffValues := RS_ZeroVector(CC,N);
			g_num := CCxy!Numerator(Rep);
			g_den := CCxy!Denominator(Rep);
			for j in [1..#Abscissas] do
				x_j := RS_PathEvaluate(Gamma,(1/2)*(Abscissas[j]+1));
				//x_j := RS_PathEvaluate(Gamma,Abscissas[j]);
				y_j := Gamma_y_Values[j];
				G_num := RS_MPolEval(g_num,x_j);
				G_den := RS_MPolEval(g_den,x_j);
				Const := RS_DerivativePathEvaluate(Gamma,(1/2)*(Abscissas[j]+1));
				//Const := RS_DerivativePathEvaluate(Gamma,Abscissas[j]);
				for Sheet in [1..N] do
					Num_Sheet := Evaluate(G_num,y_j[Sheet]);  
					Den_Sheet := Evaluate(G_den,y_j[Sheet]);
					IntegrandValue_j_Sheet := Const * (Num_Sheet/Den_Sheet);
					PathDiffValues`Entries[Sheet] +:= Weights[j]*IntegrandValue_j_Sheet;
				end for;
			end for;
			if RS_SumEntries(PathDiffValues) ge Err then
				print "RS_SumEntries(PathDiffValues):",RS_SumEntries(PathDiffValues);
				error Error("Error while integrating.");
			end if;
			Append(~PathDiffIntegrals,PathDiffValues);
		end for;
		Gamma`Integrals := PathDiffIntegrals;
	end for;


	// Arrays		
	ClosedChains := [];
	LocalMonodromy := [];
	BranchPoints := <>;

	// Construct chains from path pieces
	for l in [1..#IndexPathLists] do
		IndexList := IndexPathLists[l];
		NextPath := [];
		for Index in IndexList do
			if Sign(Index) eq 1 then
				Append(~NextPath,PathPieces[Index]);
			else
				Append(~NextPath,RS_ReversePath(PathPieces[-Index]));
			end if;
		end for;
		NextChain := RS_Chain(NextPath:Center:=DiscriminantPoints[l]);
		if NextChain`Permutation ne Id then
			Append(~ClosedChains,NextChain);
			Append(~BranchPoints,NextChain`Center);
		else
			vprint PM,2 : "Singularity detected: Discriminant point nr.",l,"is not a branch point:";
			vprint PM,2 : C_20!OrdDiscPts[l];
		end if;
	end for;

	// Get local monodromy from chains 
	LocalMonodromy := [ Ch`Permutation : Ch in ClosedChains ];
		
	// Compute permutation at infinity by relation
	InfPerm := Inverse(&*LocalMonodromy);

	// Chain around infinity
	if InfPerm ne Id then
		InfChain := (&*[ Ch : Ch in ClosedChains ])^-1;
		InfChain`Center := Infinity();
		Append(~ClosedChains,InfChain);
		Append(~LocalMonodromy,InfPerm);
		Append(~BranchPoints,Infinity());
	end if;

	vprint PM,2 : "Local monodromy:",LocalMonodromy;
	vprint PM,2 : "Number of branch points:",#BranchPoints;

	// Compute the homology
	Alpha, Cycles, BranchPlaces := RS_Tretkoff( LocalMonodromy : Genus := g );

	assert NumberOfCycles eq #Cycles;
	vprint PM,2 : "Number of cycles:",NumberOfCycles;
	vprint PM,2 : "Cycles:",Cycles;
	vprint PM,2 : "Symplectic transformation:",Alpha;


	PM_A := ZeroMatrix(CC,g,g);
	PM_B := ZeroMatrix(CC,g,g);
	
	for j in [1..g] do
		CycleIntegralValues := [];
		for k in [1..NumberOfCycles] do
			Cycle := Cycles[k];
			Value := Zero(CC);
			l := 1;
			while l lt #Cycle do
				Sheet1 := Cycle[l];
				Sheet2 := Cycle[l+2];
				PathNumber := BranchPlaces[Cycle[l+1]][1];
				Perm := LocalMonodromy[PathNumber];	
				while Sheet1 ne Sheet2 do
					Value +:= ClosedChains[PathNumber]`Integrals[j]`Entries[Sheet1];
					Sheet1 := Sheet1^Perm;
				end while;
				l +:= 2;
			end while;
			Append(~CycleIntegralValues,Value);
		end for;
		vprint PM,3 : "CycleIntegralValues:",CycleIntegralValues;	
		for m in [1..g] do
			PM_A[j][m] := &+[ Alpha[m][l] * CycleIntegralValues[l] : l in [1..NumberOfCycles] ];
			PM_B[j][m] := &+[ Alpha[m+g][l] * CycleIntegralValues[l] : l in [1..NumberOfCycles] ];
		end for;			
	end for;


	// Testing the dependence of the cycles 2g+1,...,2g+N-1
	for k in [2*g+1..NumberOfCycles] do
		CycleSum_k := &+[ Alpha[k][j] * CycleIntegralValues[j] : j in [1..NumberOfCycles] ];
		vprint PM,2 : "Cycle sum nr.",k-2*g,":",CycleSum_k;
		//if Abs(CycleSum_k) ge Err then
		if Abs(CycleSum_k) ge 10^-10 then
		//if Abs(CycleSum_k) ge 10^-(Prec+1) then
			print "Cycles summed over row nr.",k,"do not add up to zero: Computation failed.";
			print "DefiningPolynomial:",f;
			RS_SetGlobalPrecision(GlobalPrec);
			error Error("Cycles do not add up to zero: Homology test failed.");
		end if;
	end for;
	

	//vprint PM,3 : "PM_A:",PM_A;
	//vprint PM,3 : "PM_B:",PM_B;

	if Big then
		// Compute big period matrix
		PM := HorizontalJoin(PM_A,PM_B) + ZeroMatrix(C_,g,2*g);
		vprint PM,1 : "Period matrix:";
		RS_SetGlobalPrecision(GlobalPrec);
		return PM;
	else
		
		// Compute small period matrix
		PM := PM_A^(-1) * PM_B;

		// Prepare output
		MaxEntry := 0;
		for j in [1..g] do
			for k in [1..g] do
				if Abs(Re(PM[j][k])) lt RS_GetGlobalError() then
					PM[j][k] := CC_0 + Im(PM[j][k])*i;
				end if;
				if Abs(Im(PM[j][k])) lt RS_GetGlobalError() then
					PM[j][k] := Re(PM[j][k]) + CC_0*i;		
				end if;
				MaxEntry := Max(MaxEntry,Abs(PM[j][k]));
			end for;
		end for;
		PM := ChangeRing(PM,C_);


		// Testing for symmetry of the period matrix
		MaxSymDiff := 0;
		for j in [1..g] do
			for k in [j+1..g] do
				MaxSymDiff := Max(MaxSymDiff,Abs(PM[j][k] - PM[k][j]));
			end for;
		end for;
		vprint PM,1 : "Maximal symmetry deviation:",MaxSymDiff;
	
		if MaxSymDiff ge Err then
			//print "Period matrix not symmetric: Computation failed.";
			RS_SetGlobalPrecision(GlobalPrec);
			print "DefiningPolynomial:",f;
			error Error("Period matrix not symmetric: Computation failed.");
			return PM;
		end if;	
		
		// Testing positive definiteness of the imaginary part of the period matrix
		PM_Im := ZeroMatrix(RealField(Prec),g,g);
		for j in [1..g] do
			for k in [j..g] do
				PM_Im[j][k] := Real(Im(PM[j][k]));
				PM_Im[k][j] := PM_Im[j][k];
			end for;
		end for;
		assert IsPositiveDefinite(PM_Im);
		RS_SetGlobalPrecision(GlobalPrec);
		vprint PM,1 : "Period matrix:";
		return PM;
	end if;
	
end intrinsic;
