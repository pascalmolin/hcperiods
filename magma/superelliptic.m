//Import global settings
import "integration.m": RS_IntParameters;
import "analyticcontinuation.m": RS_ACRecursion;
import "pathmethods.m":RS_FindSet;
import "comparefunctions.m":RS_CompareFldComElt;

declare verbose SE,3;
declare verbose SEP,2;
declare verbose SAC,1;

C_20<i> := ComplexField(20);
C_Pi := Real(Pi(C_20));

function RS_MoebiusTransform( Points : Ht := 1000 )
	i := Parent(Points[1]).1;
	a := Random([-Ht..Ht])/Random([1..Ht]) + i*Random([-Ht..Ht])/Random([1..Ht]);
	b := Random([-Ht..Ht])/Random([1..Ht]) + i*Random([-Ht..Ht])/Random([1..Ht]);
	c := Random([-Ht..Ht])/Random([1..Ht]) + i*Random([-Ht..Ht])/Random([1..Ht]);
	d := Random([-Ht..Ht])/Random([1..Ht]) + i*Random([-Ht..Ht])/Random([1..Ht]);
	if a*d-b*c eq 0 then
		return RS_MoebiusTransform( Points );
	else
		return [ (a*x+b)/(c*x+d) : x in Points ];
	end if;
end function;

function RS_MakeCCVectors( i1, i2, Points )
	assert Type(i1) eq RngIntElt;
	a := Points[i1];
	if Type(i2) eq RngIntElt then
		b := Points[i2];
		if i1 lt i2 then
			Pts := Remove(Remove(Points,i1),i2-1);
		else
			Pts := Remove(Remove(Points,i2),i1-1);
		end if;
		CCV := [ (2*Pts[k]-b-a)/(b-a) : k in [1..#Pts] ];
	elif Type(i2) eq FldComElt then
		p := i2;
		Pts := Remove(Points,i1);
		CCV := [ (2*Pts[k]-p-a)/(p-a) : k in [1..#Pts] ];
	else
		error Error("Not supposed to happen.");
	end if;
	return CCV;
end function;
function RS_SETau3( P1,P2,P3 : Lambda := C_Pi/2 )
// Computes the maximal area of holomorphy between P1,P2 and P3
	CC := (2*P3 - P2 - P1) / (P2 - P1); // Circumcenter of the triangle defined by P1, P2, P3
	Tau_CC := Argsinh(Argtanh(CC)/Lambda);
	Tau_3 := Abs(Im(Tau_CC));
	return Tau_3;
end function;
procedure RS_SETau( Points, ~Edge )
// Computes Tau for the line segment Gamma
	for k in [1..#Points] do
		if &and[k ne Edge[2], k ne Edge[3]] then
			Edge[1] := Min(Edge[1],RS_SETau3(Points[Edge[2]],Points[Edge[3]],Points[k]));
		end if;
	end for;
end procedure;
function RS_SETau2( Points, Edge )
// Computes Tau for the line segment Gamma
	Tau := 4.0;
	for k in [1..#Points] do
		if &and[k ne Edge[2], k ne Edge[3]] then
			Tau := Min(Tau,RS_SETau3(Points[Edge[2]],Points[Edge[3]],Points[k]));
		end if;
	end for;
	return Tau;
end function;
function RS_SEInfTau( P, Points : Lambda := C_Pi/2 )
// Computes Tau for the line to infinity
	Tau := 4.0;
	for k in [1..#Points] do
		if Points[k] ne 0 then
			z := 1 - (2*P/Points[k]);
			Tau_z := Abs(Im(Argsinh(Argtanh(z)/Lambda)));
			Tau := Min(Tau,Tau_z);
		end if;
	end for;
	return Tau;
end function;
function RS_SEMST2( Points )
	d := #Points;
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
	
	// Sort by holomorphicity
	Sort(~Edges);
	
	Taken := [ 0 : j in [1..d] ];
	MST_Edges := [<Edges[1][2],Edges[1][3]>];
	Tau := Edges[1][1];
	Taken[Edges[1][2]] := 1;
	Taken[Edges[1][3]] := 1;
	while #MST_Edges lt d-1 do
		j := 2;
		while Taken[Edges[j][2]] eq Taken[Edges[j][3]] do
			j +:= 1;
		end while;
		if Taken[Edges[j][3]] eq 1 then
			Append(~MST_Edges,<Edges[j][3],Edges[j][2]>);
		else
			Append(~MST_Edges,<Edges[j][2],Edges[j][3]>);
		end if;
		Tau := Max(Tau,Edges[j][1]);
		Taken[Edges[j][2]] := 1;
		Taken[Edges[j][3]] := 1;
	end while;		

	assert #MST_Edges eq d-1;
	return MST_Edges, -Tau;
end function;
function RS_SEMST( Points )
// Computes a minimal spanning tree w.r.t holomorphicy
	d := #Points;
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
	
	// Sort by holomorphicity
	Sort(~Edges);
	Tau := -4.0;
	//print "Edges:",Edges;
	// Construct MST
	Sets := [ {v} : v in [1..d] ];
	MST_Edges := []; j := 1;
	while #MST_Edges lt d-1 do
		S1, s1 := RS_FindSet(Edges[j][2],Sets);
		S2, s2 := RS_FindSet(Edges[j][3],Sets);
		if S1 ne S2 then
			Append(~MST_Edges,<Edges[j][2],Edges[j][3]>);
			Tau := Max(Tau,Edges[j][1]);
			Sets[s1] join:= Sets[s2];
			Remove(~Sets,s2);
		end if;
		j +:= 1;
	end while;

	// Sort by increasing starting point
	Sort(~MST_Edges);

	// Order of edges in MST should suggest a path
	OrdMSTEdges := [ MST_Edges[1] ];  Taken := { MST_Edges[1][1],MST_Edges[1][2] }; Remove(~MST_Edges,1);
	while #MST_Edges ne 0 do
		for j in [1..#MST_Edges] do
			if MST_Edges[j][1] in Taken then
				Append(~OrdMSTEdges,MST_Edges[j]);
				Taken join:= {MST_Edges[j][2]};
				Remove(~MST_Edges,j);				
				break;
			elif MST_Edges[j][2] in Taken then
				Append(~OrdMSTEdges,<MST_Edges[j][2],MST_Edges[j][1]>);
				Taken join:= {MST_Edges[j][1]};
				Remove(~MST_Edges,j);				
				break;
			end if;
		end for;
	end while;

	assert #OrdMSTEdges eq d-1;
	return OrdMSTEdges, -Tau;
end function;
function RS_SEMST2( Points : ED := false)
	d := #Points;
	Edges := [];
	for k in [1..d] do
		for l in [k+1..d] do
			if ED then
				Edge := <-Abs(Points[k]-Points[l]),k,l>;
			else
				Edge := <4.0,k,l>;
				RS_SETau( Points, ~Edge );
			end if;
			Append(~Edges,Edge);
		end for;
	end for;
	
	// Sort by holomorphicity
	Sort(~Edges);
	Taken := [ 0 : j in [1..d] ];
	k := 1; n := #Edges;
	if ED then
		Tau := RS_SETau2(Points,Edges[n]);
	else
		Tau := Edges[n][1];
	end if;
	MST_Edges := [ <Edges[n][2],Edges[n][3]> ];
	Taken[Edges[n][2]] := 1;
	Taken[Edges[n][3]] := 1;
	while k lt d-1 do
		l := n-1;
		while Taken[Edges[l][2]] eq Taken[Edges[l][3]] do
			l -:= 1;
		end while;
		if Taken[Edges[l][3]] eq 1 then
			// Flip edge
			Append(~MST_Edges,<Edges[l][3],Edges[l][2]>);
		else
			Append(~MST_Edges,<Edges[l][2],Edges[l][3]>);	
		end if;
		Taken[Edges[l][2]] := 1;
		Taken[Edges[l][3]] := 1;
		if ED then
			Tau := Min(Tau,RS_SETau2(Points,Edges[l]));
		else
			Tau := Min(Tau,Edges[l][1]);
		end if;
		k +:= 1;
	end while;
	assert #MST_Edges eq d-1;
	return MST_Edges, Tau;
end function;


function RS_ImSgn( z )
	if IsReal(z) then
		return 1;
	else
		return Sign(Im(z));
	end if;
end function;


function RS_NRootAC(x,CCV,Zetas,N:long:=false)
// Analytic continuation of y = f(x)^(1/N) on one sheet, avoiding branch cuts
	//long := true;
	/*if #CCV eq 1 or long then
		prod := 1;
		for k in [1..#CCV] do
			if Real(CCV[k]) gt 0 then
				prod *:= (CCV[k] - x)^(1/N);
			else
				prod *:= (x - CCV[k])^(1/N);
			end if;
		end for;
		return prod;
	end if;*/
	WindingNr := 0;
	if Re(CCV[1]) le 0 then
		z := x - CCV[1];
	else
		z := CCV[1] - x;
	end if;
	isgnz := RS_ImSgn(z);
	for k in [2..#CCV] do
		if Re(CCV[k]) le 0 then
			z_ := x - CCV[k];
		else
			z_ := CCV[k] - x;
		end if;
		isgnz_ := RS_ImSgn(z_);
		z *:= z_;
		if isgnz eq isgnz_ then
			isgnz := RS_ImSgn(z);
			if isgnz ne isgnz_ then
				if isgnz_ gt 0 then
					WindingNr +:= 2;
				else
					WindingNr -:= 2;
				end if;
			end if;
		else
			isgnz := RS_ImSgn(z);
		end if;
	end for;
	
	if Im(z) eq 0 and Re(z) lt 0 then
		z := -z;
		WindingNr +:= 1;	
	end if;
	
	if N eq 2 then
		z := Zetas[WindingNr mod (2*N) + 1]  * Sqrt(z);
	else
		z := Zetas[WindingNr mod (2*N) + 1]  * Exp((1/N)*Log(z));
	end if;
	return z;
end function;



function RS_SEIntegrate( Edge,Points,DFF,Abscissas,Weights,StepLength,d,N,Zetas,lc_m )

	CC<i> := Parent(Points[1]); CC_0 := Zero(CC); g := &+[ #DFF_i : DFF_i in DFF ];

	a := Points[Edge[1]];
	b := Points[Edge[2]];
	
	// Make vector of centers of circumcircles
	CCV := RS_MakeCCVectors(Edge[1],Edge[2],Points);

	// Factors due to change of variable
	Fact1 := (b-a)/2;
	Fact2 := (b+a)/(b-a);
	UP := #[ z : z in CCV | Re(z) gt 0 ];

	vprint SE,2 : "################### Next Integral ###################";
	Integrals0 := [ CC_0 : j in [1..g] ];

	Integrals := []; // g x N array of integrals
	ElementaryIntegrals := []; // Needed for Abel-Jacobi

	// Initiate on x = 0, dx = 1
	z := 1/RS_NRootAC(0,CCV,Zetas,N);
	Enum := Fact2^DFF[1][1][1];
	Denom := z^DFF[1][1][2];
	Denomi := Denom;
	Ct := 1;
	Integrals0[Ct] +:= (Enum * Denomi);
	for k in [2..#DFF[1]] do
		Denomi *:= z;
		Ct +:= 1;
		Integrals0[Ct] +:= (Enum * Denomi);
	end for;
	for j in [2..#DFF] do
		Enum *:= Fact2;
		Diff := DFF[j][1][2] - DFF[j-1][1][2];
		for l in [1..Diff] do
			Denom *:= z;
		end for;
		Denomi := Denom;
		Ct +:= 1;
		Integrals0[Ct] +:= (Enum * Denomi);
		for k in [2..#DFF[j]] do		
			Denomi *:= z;
			Ct +:= 1;
			Integrals0[Ct] +:= (Enum * Denomi);
		end for;
	end for;
	
	for t in [1..#Abscissas] do
		x := Abscissas[t];
		// Analytic continuation
		z1 := 1/RS_NRootAC(x,CCV,Zetas,N);
		z2 := 1/RS_NRootAC(-x,CCV,Zetas,N);
		x1 := x + Fact2;
		x2 := -x + Fact2;
		Enum1 := x1^DFF[1][1][1];
		Enum2 := x2^DFF[1][1][1];
		Denom1 := z1^DFF[1][1][2];
		Denom2 := z2^DFF[1][1][2];
		Denomi1 := Denom1;
		Denomi2 := Denom2;
	
		// Weights
		w2_2N_inv := Exp( (2/N) * -Log(Weights[t][2]));
		dx := Weights[t][1] * Weights[t][2]^2 * w2_2N_inv^DFF[1][1][2];
		dx_ := dx;

		Ct := 1;
		Integrals0[Ct] +:= (Enum1 * Denomi1 + Enum2 * Denomi2) * dx_;
		for k in [2..#DFF[1]] do
			Denomi1 *:= z1;
			Denomi2 *:= z2;
			Ct +:= 1;
			dx_ *:= w2_2N_inv;
			Integrals0[Ct] +:= (Enum1 * Denomi1 + Enum2 * Denomi2) * dx_;
		end for;

		for j in [2..#DFF] do
			Enum1 *:= x1;
			Enum2 *:= x2;
			Diff := DFF[j][1][2] - DFF[j-1][1][2];
			for l in [1..Diff] do
				Denom1 *:= z1;
				Denom2 *:= z2;
				dx *:= w2_2N_inv;
			end for;
			Denomi1 := Denom1;
			Denomi2 := Denom2;
			dx_ := dx;
			Ct +:= 1;
			Integrals0[Ct] +:= (Enum1 * Denomi1 + Enum2 * Denomi2) * dx_;
			for k in [2..#DFF[j]] do
				Denomi1 *:= z1;
				Denomi2 *:= z2;
				Ct +:= 1;
				dx_ *:= w2_2N_inv;
				Integrals0[Ct] +:= (Enum1 * Denomi1 + Enum2 * Denomi2) * dx_;
			end for;
		end for;
	end for;

	
	Factor1 := StepLength * Fact1^(DFF[1][1][1]+1); // Lambda * h * ((b-a)/2)^(i+1);
	Fact1_dN_inv := Exp( (d/N) * -Log(Fact1)); // (b-a)^-(d/N)
	Factor2 := Fact1_dN_inv^DFF[1][1][2]; // (b-a)^-(j*d/N)
	Factor21 := Factor2;
	Ct := 1;
	Pow := -((UP + 1) mod 2)*DFF[1][1][2] mod (2*N) + 1;
	Integrals0[Ct] *:= lc_m^DFF[1][1][2] * Zetas[Pow] * Factor1 * Factor21;
	Append(~ElementaryIntegrals,Integrals0[Ct]);
	Pow := -2*DFF[1][1][2] mod (2*N) + 1;
	Integrals0[Ct] *:= (1 - Zetas[Pow]);
	NextIntegrals := [ Integrals0[Ct] ];
	for l in [1..N-2] do
		Pow := -2*l*DFF[1][1][2] mod (2*N) + 1;
		Append(~NextIntegrals,Zetas[Pow]*Integrals0[Ct]);
	end for; 
	Append(~Integrals,NextIntegrals);
	for k in [2..#DFF[1]] do
		Factor21 *:= Fact1_dN_inv;
		Ct +:= 1;
		Pow := -((UP + 1) mod 2)*DFF[1][k][2] mod (2*N) + 1;
		Integrals0[Ct] *:= lc_m^DFF[1][k][2] * Zetas[Pow] * Factor1 * Factor21;
		Append(~ElementaryIntegrals,Integrals0[Ct]);
		Pow := -2*DFF[1][k][2] mod (2*N) + 1;
		Integrals0[Ct] *:= (1 - Zetas[Pow]);
		NextIntegrals := [ Integrals0[Ct] ];
		for l in [1..N-2] do
			Pow := -2*l*DFF[1][k][2] mod (2*N) + 1;
			Append(~NextIntegrals,Zetas[Pow]*Integrals0[Ct]);
		end for; 
		Append(~Integrals,NextIntegrals);
	end for;

	for j in [2..#DFF] do
		Factor1 *:= Fact1;
		Diff := DFF[j][1][2] - DFF[j-1][1][2];
		for l in [1..Diff] do
			Factor2 *:= Fact1_dN_inv;
		end for;
		Factor21 := Factor2;
		Ct +:= 1;
		Pow := -((UP + 1) mod 2)*DFF[j][1][2] mod (2*N) + 1;
		Integrals0[Ct] *:= lc_m^DFF[j][1][2] * Zetas[Pow] * Factor1 * Factor21;
		Append(~ElementaryIntegrals,Integrals0[Ct]);
		Pow := -2*DFF[j][1][2] mod (2*N) + 1;
		Integrals0[Ct] *:= (1 - Zetas[Pow]);
		NextIntegrals := [ Integrals0[Ct] ];
		for l in [1..N-2] do
			Pow := -2*l*DFF[j][1][2] mod (2*N) + 1;
			Append(~NextIntegrals,Zetas[Pow]*Integrals0[Ct]);
		end for; 
		Append(~Integrals,NextIntegrals);
		for k in [2..#DFF[j]] do
			Factor21 *:= Fact1_dN_inv;
			Ct +:= 1;
			Pow := -((UP + 1) mod 2)*DFF[j][k][2] mod (2*N) + 1;
			Integrals0[Ct] *:= lc_m^DFF[j][k][2] * Zetas[Pow] * Factor1 * Factor21;
			Append(~ElementaryIntegrals,Integrals0[Ct]);
			Pow := -2*DFF[j][k][2] mod (2*N) + 1;
			Integrals0[Ct] *:= (1 - Zetas[Pow]);
			NextIntegrals := [ Integrals0[Ct] ];
			for l in [1..N-2] do
				Pow := -2*l*DFF[j][k][2] mod (2*N) + 1;
				Append(~NextIntegrals,Zetas[Pow]*Integrals0[Ct]);
			end for; 
			Append(~Integrals,NextIntegrals);
		end for;
	end for;
	return Integrals, ElementaryIntegrals;
end function;


function RS_SEIntegrate2( Edge,Points,DFF,Abscissas,Weights,StepLength,d,N,Zetas,lc_m )

	CC<i> := Parent(Points[1]); CC_0 := Zero(CC); g := &+[ #DFF_i : DFF_i in DFF ];
	a := Points[Edge[1]];
	b := Points[Edge[2]];
	// Make vector of centers of circumcircles
	CCV := RS_MakeCCVectors(Edge[1],Edge[2],Points);

	// Factors due to change of variable
	Fact1 := (b-a)/2;
	Fact2 := (b+a)/(b-a);
	UP := #[ z : z in CCV | Re(z) gt 0 ];
	vprint SE,2 : "################### Next Integral ###################";
	Integrals0 := [ CC_0 : j in [1..g] ];

	Integrals := []; // g x N array of integrals
	ElementaryIntegrals := []; // Needed for Abel-Jacobi
	
	DFF := &cat[ DFF_i : DFF_i in DFF ];
	g := #DFF;

	// Initiate on x = 0, dx = 1
	z := RS_NRootAC(0,CCV,Zetas,N);
	for j in [1..g] do
		w := DFF[j];
		Denom := z^w[2];
		Integrals0[j] +:= (Fact2^w[1]/Denom);
	end for;

	for t in [1..#Abscissas] do
		x := Abscissas[t];
		
		// Analytic continuation
		z1 := RS_NRootAC(x,CCV,Zetas,N);
		z2 := RS_NRootAC(-x,CCV,Zetas,N);

		Enum1 := (x + Fact2);
		Enum2 := (-x + Fact2);

		for j in [1..g] do
			w := DFF[j];
			Denom1 := z1^w[2];
			Denom2 := z2^w[2];
			dx := Weights[t][1] * Weights[t][2]^(2-(2*w[2]/N));
			Integrals0[j] +:= (Enum1^w[1]/Denom1 + Enum2^w[1]/Denom2) * dx;
		end for;
	end for;
	
	for j in [1..g] do
		w := DFF[j];
		Factor := Fact1^((w[1]+1)-(d*w[2]/N));
		Pow := -((UP + 1) mod 2)*w[2] mod (2*N) + 1;
		Integrals0[j] *:= lc_m^w[2] * Zetas[Pow] * StepLength * Factor;
		Append(~ElementaryIntegrals, Integrals0[j]);
		Integrals0[j] *:= ( 1 - Zetas[3]^-w[2]);
		Append(~Integrals,[ Zetas[3]^(-l*w[2]) * Integrals0[j] : l in [0..N-2]]);
	end for;
	return Integrals, ElementaryIntegrals;
end function;


procedure RS_PolyShift( ~Vec, Const )
	g := #Vec;
	for k in [1..g] do;
		l := g;
		while l ge k do
			Vec[l] +:= Vec[l-1] * Const;
			l -:= 1;
		end while;
	end for;
end procedure;
function RS_SEIntegrate3( Edge,Points,DFF,Abscissas,Weights,StepLength,d,N,Zetas,lc_m )

	CC<i> := Parent(Points[1]); CC_0 := Zero(CC); g := &+[ #DFF_i : DFF_i in DFF ];
	a := Points[Edge[1]];
	b := Points[Edge[2]];
	// Make vector of centers of circumcircles
	CCV := RS_MakeCCVectors(Edge[1],Edge[2],Points);
	//UP := (#[ z : z in CCV | Re(z) gt 0 ] + 1) mod 2;
	UP := #[ z : z in CCV | Re(z) gt 0 ];

	// ((b-a)/2)^i, i = 1,..,d
	z := (b-a)/2; Fact1 := [ z ];
	for j in [2..d] do
		Append(~Fact1,z*Fact1[j-1]);
	end for;
	// (2/(b-a))^(dj/N), j = 1,..,N-1
	z := Exp( (d/N) * -Log(z) ); Fact2 := [ z ];
	for j in [2..N-1] do
		Append(~Fact2,z*Fact2[j-1]);
	end for;
	//
	z := (b+a)/(b-a); Fact3 := [ 1,z ];
	for j in [3..d-1] do
		Append(~Fact3,z*Fact3[j-1]);
	end for;
	vprint SE,2 : "################### Next Integral ###################";
	
	g := #DFF;
	ShiftedIntegrals := [ CC_0 : j in [1..g] ];


	// Initiate on x = 0, dx = 1
	y1 := 1/RS_NRootAC(0,CCV,Zetas,N); // 1/y(0)
	wy1 := y1;
	iy := 1;
	for j in [1..g] do
		w := DFF[j];
		while iy lt w[2] do
			wy1 *:= y1;
			iy +:= 1;
		end while;
		//ShiftedIntegrals[j] +:= wy1*z^w[1];
		if w[1] eq 0 then
			ShiftedIntegrals[j] +:= wy1;
		end if;
	end for;

	for t in [1..#Abscissas] do
		x := Abscissas[t];
		y1 := 1/RS_NRootAC(x,CCV,Zetas,N); // 1/y(x)
		y2 := 1/RS_NRootAC(-x,CCV,Zetas,N); // 1/y(-x)
		w2N := Exp( -(2/N) * Log(Weights[t][2])); // cosh^-(2/N)
		y1 *:= w2N;
		y2 *:= w2N; 
		w12 := Weights[t][1] * Weights[t][2]^2;
		wy1 := y1 * w12;
		wy2 := y2 * w12;
		wyx1 := wy1;
		wyx2 := wy2;
		ix := 0; iy := 1;
		for j in [1..g] do
			w := DFF[j];
			while iy lt w[2] do
				wy1 *:= y1;
				wy2 *:= y2;
				wyx1 := wy1;
				wyx2 := wy2;
				iy +:= 1;
				ix := 0;
			end while;
			while ix lt w[1] do
				wyx1 *:= x;
				wyx2 *:= -x;
				ix +:= 1;
			end while;
			ShiftedIntegrals[j] +:= wyx1 + wyx2;
		end for;
	end for;

	Integrals0 := [ CC_0 : j in [1..g] ];
	ElementaryIntegrals := []; // Needed for Abel-Jacobi
	Integrals := []; // g x N array of integrals

	for j in [1..g] do
		w := DFF[j];
		for k in [0..w[1]] do
			print "Binomial(w[1],k) :",Binomial(w[1],k) ;
			Integrals0[j] +:= Binomial(w[1],k) * Fact3[k+1] * ShiftedIntegrals[j-k];
		end for;
		Pow := -((UP + 1) mod 2)*DFF[j][2] mod (2*N) + 1;
		Integrals0[j] *:= lc_m^w[2] * Zetas[Pow] * StepLength * Fact1[w[1]+1] * Fact2[w[2]];
		Append(~ElementaryIntegrals, Integrals0[j]);
		Pow := -2*DFF[j][2] mod (2*N) + 1;
		Integrals0[j] *:= (1 - Zetas[Pow]);
		NextIntegrals := [ Integrals0[j] ];
		for l in [1..N-2] do
			Pow := -2*l*DFF[j][2] mod (2*N) + 1;
			Append(~NextIntegrals,Zetas[Pow]*Integrals0[j]);
		end for; 
		Append(~Integrals,NextIntegrals);
	end for;
	return Integrals, ElementaryIntegrals;
end function;


intrinsic RS_SEPM( f::RngMPolElt : Prec := -1, Small := true ) -> Mtrx
{ Computes period matrices associated to the superelliptic curve defined by f to precision Prec }
	
	vprint SE,1 : "Defining polynomial:",f;

	// Precision
	if Prec lt 20 then
		Prec := 20;
	end if;
	vprint SE,1 : "Prescribed precision:",Prec;
	
	// Error
	Err := 10^-(Prec+1);

	// Degrees and genus
	m := Degree(f,2); n := Degree(f,1); g := Round((1/2) * ((m-1)*(n-1) - Gcd(m,n) + 1));
	assert &and[m ge 2, n ge 3];
	vprint SE,1 : "m:",m;
	vprint SE,1 : "n:",n;
	vprint SE,1 : "Genus:",g;


	// Complex fields and constants
	C_<i> := ComplexField(Prec);
	CC<i> := ComplexField(Prec+10); CCx<x> := PolynomialRing(CC);
	CC_0 := Zero(CC); CC_Pi := Pi(CC);
	
	vprint SEP,1 : "Precision CC_0:",Precision(CC_0);
	vprint SEP,1 : "Precision Pi:",Precision(CC_Pi);

	// Root of unity powers
	Zetas := [ RS_ReducePrecision(Exp(k*CC_Pi*i/m),5) : k in [0..2*m-1] ];
	Zeta := Zetas[3];
	vprint SE,1 : "Zeta:",C_20!Zeta;
	vprint SEP,1 : "Precision Zeta:",Precision(Zeta);

	// Branch Points
	f_x := CCx!UnivariatePolynomial(Evaluate(f,[Parent(f).1,0]));
	Points :=  RS_Roots(f_x);

	// Leading coefficient
	lc_x := -LeadingCoefficient(f_x); 
	vprint SE,2 : "Leading coefficient:",C_20!lc_x;	
	lc_m := 1/lc_x^(1/m);

	// Low precision branch points
	LowPrecPoints := [ C_20!Pt : Pt in Points ];

	vprint SE,2 : "Discriminant points:",LowPrecPoints;
	vprint SEP,1 : "Precision of discriminant points:",Precision(Points[1]);;

	// Maximal spanning tree
	vprint SE,1 : "Constructing spanning tree...";
	t := Cputime();
	if n ge 15 then
		vprint SE,1 : "using euclidean weights...";
		MST_Edges, MST_Tau := RS_SEMST2(LowPrecPoints:ED:=true);
	else
		vprint SE,1 : "using holomorphic weights...";
		MST_Edges, MST_Tau := RS_SEMST2(LowPrecPoints:ED:=false);
	end if;
	Cputime(t);
	/*
	vprint SE,1 : "Improving Tau...",C_20!MST_Tau;
	t := Cputime();
	while MST_Tau lt 0.2 do
		vprint SE,1 : "Tau too bad: applying moebius transform...";
		Points := RS_MoebiusTransform(Points); 
		LowPrecPoints := [ C_20!Pt : Pt in Points ];
		MST_Edges, MST_Tau := RS_SEMST2(LowPrecPoints:ED);
		vprint SE,1 : "MST_Tau:",C_20!MST_Tau;
	end while;
        Cputime(t);
	*/
	vprint SE,2 : "MST_Edges:",MST_Edges;
	vprint SE,2 : "#MST_Edges:",#MST_Edges;
	vprint SE,1 : "MST_Tau:",C_20!MST_Tau;

	// Holomorphic differentials
	vprint SE,1 : "Computing holomorphic differentials...";
	DFF := RS_SEDifferentials(m,n:Int3:=false);
	vprint SE,2 : "Holomorphic differentials:",DFF;

	// Integration parameters
	vprint SE,1 : "Computing Integration parameters...";
	t := Cputime();
	Abscissas, Weights, StepLength := RS_SEDEIntegrationParameters(LowPrecPoints,MST_Edges,MST_Tau,m,CC);
	Cputime(t);

	// Computing periods
	vprint SE,1 : "Integrating...";
	t := Cputime();
	Integrals := [];
	for k in [1..n-1] do
		vprint SE,2 : "Integrating edge nr.",k;
		I := RS_SEIntegrate( MST_Edges[k],Points,DFF,Abscissas,Weights,StepLength,n,m,Zetas,lc_m );
		//I := RS_SEIntegrate3( MST_Edges[k],Points,DFF,Abscissas,Weights,StepLength,n,m,Zetas,lc_m );
		Append(~Integrals,I);
	end for;
	Cputime(t);

	// Period matrix
	PM_ := ZeroMatrix(CC,g,(m-1)*(n-1));
	for k in [1..n-1] do
		for l in [1..m-1] do
			Ind := (m-1)*(k-1) + l;
			for j in [1..g] do
				PM_[j][Ind] := Integrals[k][j][l];
			end for;
		end for;
	end for;
	vprint SE,3: "Integrals:",PM_;


	// Intersection matrix
	vprint SE,1 : "Computing skd-matrix...";
	t := Cputime();
	spsm_Matrix := [ [] : j in [1..n-1] ];
	for j in [1..n-1] do
		spsm_Matrix[j][j] := <1,m-1>;
		for l in [j+1..n-1] do
			spsm := RS_SEIntersection(MST_Edges[j],MST_Edges[l],Points,m,Zetas);
			vprint SE,2: "<sp,sm> pair of edges",MST_Edges[j],"and",MST_Edges[l],":",<spsm[1] mod m,spsm[2] mod m>;
			spsm_Matrix[j][l] := <spsm[1] mod m,spsm[2] mod m>;
		end for;
	end for;
	Cputime(t);
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
	RealIntMatrix := BlockMatrix(n-1,n-1,Blocks);
	Cputime(t);
	vprint SE,2: "(Real) Intersection matrix:",RealIntMatrix;
	vprint SE,2: "Rank(IntersectionMatrix):",Rank(RealIntMatrix);
	assert Rank(RealIntMatrix) eq 2*g;

	// Symplectic reduction of intersection matrix
	vprint SE,1 : "Performing symplectic reduction...";
	t := Cputime();
	CF, ST := RS_SymplecticBasis(RealIntMatrix);
	Cputime(t);
	vprint SE,3: "ST:",ST;
	vprint SE,3: "CF:",CF;
	ST_CC := ChangeRing(Transpose(ST),CC);


	vprint SE,1 : "Matrix multiplication 1...";
	// Can be avoided?
	t := Cputime();
	PMAPMB := PM_ * ST_CC;
	Cputime(t);
	
	
	PM_A := ColumnSubmatrixRange(PMAPMB,1,g);
	PM_B := ColumnSubmatrixRange(PMAPMB,g+1,2*g);
	vprint SE,3 : "PM_A:",PM_A;
	vprint SE,3 : "PM_B:",PM_B;
	vprint SE,3 : "Dependent columns:",ColumnSubmatrixRange(PMAPMB,2*g+1,Nrows(ST_CC));
	

	if Small eq false then
		// Compute big period matrix
		PM := HorizontalJoin(PM_A,PM_B) + ZeroMatrix(C_,g,2*g);
		vprint SE,1 : "Period matrix:";
		return PM;
	end if;

	vprint SE,1 : "Matrix inversion...";
	t := Cputime();
	PM_AInv := PM_A^(-1);
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
			if Abs(Re(PM[j][k])) lt Err then
				PM[j][k] := Zero(C_) + Im(PM[j][k])*i;
			end if;
			if Abs(Im(PM[j][k])) lt Err then
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
	
	if MaxSymDiff gt 10^-15 then
	//if MaxSymDiff gt Err then
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

	return PM;

end intrinsic;


intrinsic RS_SEDifferentials( N::RngIntElt, d::RngIntElt : Int3 := true ) -> SeqEnum[Tup]
{ Returns a basis of holomorphic differentials for superelliptic curves }
	DFF := [];
	delta := Gcd(N,d);
	if Int3 then
	    	for j in [1..N-1] do
			for i in [1..d-1] do
				if N*i le j*d-delta then
					Append(~DFF,<i-1,j>);
				end if;
			end for;
		end for;
		assert 2*#DFF eq (N-1)*(d-1)-delta+1;
	else
		for i in [0..d-2] do
			DFF_i := [];
			for j in [1..N-1] do
				if (-N*(i+1) + j*d - delta ge 0) then
					Append(~DFF_i,<i,j>);
				end if;
			end for;
			if #DFF_i ne 0 then 
				Append(~DFF,DFF_i);
			end if;
		end for;
		assert 2*&+[ #w : w in DFF] eq (N-1)*(d-1)-delta+1;
	end if;
	return DFF;
end intrinsic;


function RS_SEDEIntPoints( h, NrPoints, Lambda )
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
	/*
	Abscissas := Reverse([-zk : zk in Abscissas_]) cat [0] cat Abscissas_;
	Weights := Reverse(Weights_) cat [<1,1>] cat Weights_;
  	return Abscissas, Weights, h*Lambda;
	*/
	return Abscissas_, Weights_, h*Lambda;
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

function RS_M1Bound(CCV,N)
// Compute the bound M1 on the integrand
	B1 := 1/(&*[ RS_Dist_1(P) : P in CCV ])^(1/N);
	//print "B1:",B1;
	if B1 gt 1 then
		return B1^(N-1);
	else
		return B1;
	end if;
end function;

/*
function RS_tOpt(D,M1)
// Some value
	C<i> := ComplexField(20); C_Pi := Pi(C); Lambda := C_Pi/2;
	Alpha := (N-1)/N;
	return Argsinh( D+ Log((4^(1-Alpha) * M_Tau)/(2-2*Alpha))/(Lambda*(2-2*Alpha)));
end function;
*/

function RS_Dist_thsh( P, Tau : Lambda := C_Pi/2 )
// Strange distance function?!?
	P := Abs(Re(P)) + i*Abs(Im(P));
	x0 := 0; x1 := Argcosh(C_Pi/(2*Lambda*Sin(Tau))); // s.t. \Lambda Cosh(x)Sin(\Tau)= Pi/2;
	Phi := function(x)
		return Tanh(Lambda*Sinh(x+i*Tau));
	end function;
	ArgPhiP := function(x)
		return Arg(Cosh(x+i*Tau))-2*Arg(Cosh(Lambda*Sinh(x+i*Tau)));
	end function;
	ArgPhi := function(x)
		return Arg(P - Phi(x));
	end function;
	// x := Solve(x := x0, x1, ArgPhi(x)-ArgPhiP(x)-C_Pi/2);
	x := x0;
	Val := ArgPhi(x)-ArgPhiP(x)-C_Pi/2;
	//print "x0:",x;
	//print "x1:",x1;
	//print "Val:",Val;
	n := 5;
	for t in [0..n] do
		x_new := (t/n)*x1 + (1-(t/n))*x0;
		Val_new := ArgPhi(x_new)-ArgPhiP(x_new)-C_Pi/2;
		if Abs(Val_new) lt Abs(Val) then
			x := x_new;
			Val := Val_new;
		end if;
	end for;
	//print "x_new:",x_new;
	//print "Val_new:",Val_new;
	return Abs(P-Phi(x));
end function;


function RS_M2Bound(CCV,Tau,N)
// Compute the bound M2 on the integrand
	B2 := 1/(&*[ RS_Dist_thsh(P,Tau) : P in CCV ])^(1/N);
	//print "B2:",B2;
	if B2 gt 1 then
		return B2^(N-1);
	else
		return B2;
	end if;
end function;


function RS_PhiBound(Tau,N:Lambda := C_Pi/2 )
// Upper bound for transformation phi
	Alpha := 1/N;
	//Alpha := (N-1)/N;
	Y_Tau := Sqrt(C_Pi/2*Lambda*Sin(Tau));
	//Y_Tau := (Pi/2) * Sin(Tau)^(7/8);
	X_Tau := Sqrt(Y_Tau^2 - Lambda^2 * Sin(Tau)^2) / Tan(Tau);
	return ((1+2*Alpha)/(2*Alpha*Cos(Y_Tau)^(2*Alpha))+1/(2*Alpha*Sinh(X_Tau)^(2*Alpha)));
	// something different...
end function;

/*
function RS_hOpt(D,Tau,M2,N,Lambda)
	I := RS_PhiBound(Tau,N);
	return (2*Pi*Tau)/Log(1+2*M2*I*Exp(D));
end function;
*/

intrinsic RS_SEDEIntegrationParameters( Points::SeqEnum[FldComElt], MST::SeqEnum[Tup], Tau::FldReElt, N::RngIntElt, CC::FldCom ) -> SeqEnum[FldReElt], SeqEnum[FldReElt], FldReElt
{ Computes integration parameters for double exponential integration on the interval [-1,1] }
	assert N ge 2;
	// Pi
	CC_Pi := Real(Pi(CC));

	// Compute D
	D := Precision(CC) * Log(CC!10);
	vprint SE,1 : "D:",Ceiling(D);

	// Check Tau
	if Tau ge CC_Pi/2 then
		error Error("Not supposed to happen.");
	end if;
	assert &and[Tau lt CC_Pi/2, Tau gt 0];
	vprint Int,1 : "Tau:",C_20!Tau;

	// Make it precise
	Tau *:= (19/20);
	vprint SE,1 : "Tau for integration:",Tau;	

	// Parameter Lambda = Pi/2 for DE integration
	Lambda := CC_Pi/2;	

	// Compute bounds M1,M2
	M1 := 0; M2 := 0;
	if #MST gt 0 then
		for Edge in MST do
			// Make vector of centers of circumcircles
			CCV := RS_MakeCCVectors(Edge[1],Edge[2],Points);
			M1 := Max(M1,RS_M1Bound(CCV,N));
			M2 := Max(M2,RS_M2Bound(CCV,Tau,N));
		end for;
	else
		M1 := 10000;
		M2 := 10000;
	end if;
	M1 := Ceiling(M1); M2 := Ceiling(M2);
	vprint SE,1 : "M1:",M1;
	vprint SE,1 : "M2:",M2;
	
	// New parameters
	Alpha := 1/N;
	Tau := CC!Tau; D := CC!D;
	X_Tau := Cos(Tau) * Sqrt( (CC_Pi / (2*Lambda*Sin(Tau))) - 1 );
	B_TA := (2/Cos(Tau)) * ( (X_Tau/2) * (1/(Cos(Lambda*Sin(Tau))^(2*Alpha)) + (1/X_Tau^(2*Alpha))) ) + (1/(2*Alpha*Sinh(X_Tau)^(2*Alpha)));
	h := Real( 2 * CC_Pi * Tau /  ( D+Log(2*M2 * B_TA + 1 )));
	n := Ceiling(Argsinh((D+ Log((2^(2*Alpha+1) *M1)/Alpha )) / ( 2*Alpha*Lambda ))/h);
	
	vprint SE,1 : "h:",C_20!h;
	vprint SE,1 : "n:",n;

	// Compute integration points
	Abscissas, Weights, StepLength := RS_SEDEIntPoints(h,n,Lambda);
	
	return Abscissas, Weights, StepLength;

end intrinsic;


intrinsic RS_SEIntegration( f::RngMPolElt, Gamma::RSPath, Points::SeqEnum[FldComElt], DFF::SeqEnum[Tup] ) -> SeqEnum[FldComElt]
{ Integration on superelliptic curves }
	
	"##################### RS_SEIntegration ###################";
	// Variables
	d := Degree(f,1); assert d eq #Points;
	N := Degree(f,2); g := #DFF;
	assert #Points ge 2;

	a := Gamma`StartPt;
	b := Gamma`EndPt;

	k := Position(Points,Points[1]);
	l := Position(Points,Points[2]);
	assert k ne l;
	
	if k eq 0 then
		PP := Remove(Points,l);
	elif k lt l then
		PP := Remove(Remove(Points,l),k);
	else
		PP := Remove(Remove(Points,k),l);
	end if;

	CCV := [ (2*PP[k]-b-a)/(b-a) : k in [1..#PP] ]; 
	if #PP eq d-2 then
		cc := 0;
	else
		cc := 1;
	end if;

	Fact1 := (b-a)/2; print "Fact1:",Fact1;
	Fact2 := (b+a)/(b-a);

	RS_SetGlobalZeta(N);
	C<i> := RS_GetGlobalComplexField_Max(); Pi := RS_GetGlobalPi(); C_1 := One(C); CC_0 := Zero(C);
	Zeta := RS_GetGlobalZeta(); SqrtZeta := Sqrt(Zeta);

	// Leading coefficient
	f_x := UnivariatePolynomial(Evaluate(f,[Parent(f).1,0]));
	LC_x := 1/C!(-LeadingCoefficient(f_x))^(1/N);;
	print "LC_x:",LC_x;


	print "Gamma:",Gamma;
	print "PP:",PP;
	
	RS_TauPath(PP,Gamma);
	print "Gamma`Tau:",Gamma`Tau;
	
	Abscissas, Weights, StepLength := RS_SEDEIntegrationParameters(Points,[],Gamma`Tau,N,C);
	
	// Factors due to change of variable
	Fact1 := (b-a)/2;
	Fact2 := (b+a)/(b-a);

	Integrals0 := [ CC_0 : j in [1..g] ];
	// Initiate on x = 0, dx = 1
	z := RS_NRootAC(0,CCV,Zeta,N);
	for j in [1..g] do
		w := DFF[j];
		Denom := z^w[2];
		Integrals0[j] +:= (Fact2^w[1]/Denom);
	end for;
	
	for t in [1..#Abscissas] do
		x := Abscissas[t];
		
		// Analytic continuation
		z1 := RS_NRootAC(x,CCV,Zeta,N);
		z2 := RS_NRootAC(-x,CCV,Zeta,N);

		Enum1 := (x + Fact2);
		Enum2 := (-x + Fact2);

		for j in [1..g] do
			w := DFF[j];
			Denom1 := z1^w[2];
			Denom2 := z2^w[2];
			dx := Weights[t][1] * Weights[t][2]^(2-(2*w[2]/N));
			Integrals0[j] +:= (Enum1^w[1]/Denom1 + Enum2^w[1]/Denom2) * dx;
		end for;
	end for;

	Integrals := []; // g x N array of integrals
	
	for j in [1..g] do
		w := DFF[j];
		Factor := Fact1^((w[1]+1)-(d*w[2]/N));
		Integrals0[j] *:= (SqrtZeta^(-w[2]) - SqrtZeta^(-w[2])) * StepLength * Factor;
		//Integrals0[j] *:= (1 - Zeta^-w[2]) * StepLength * Factor;
		Append(~Integrals,[ Zeta^(l*w[2]) * LC_x^j * Integrals0[j] : l in [0..N-1]]);
	end for;

	return Integrals;

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
	DFF := RS_SEDifferentials(f);
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
