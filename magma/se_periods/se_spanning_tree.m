/*******************************************************************************

 Copyright (C) 2017 
 Code adapted from Pascal Molin,
 written by Christian Neurohr

*******************************************************************************/

// Import functions
import "se_help_funcs.m": MakeCCVector, SortByRealPart, IntGroups;
import "se_de_int.m": Bound_M1, Bound_M2;
import "se_gj_int.m": DistanceEllipse;


// Complex field
C<I> := ComplexField(20);
RPI := Real(Pi(C));

// Compare 3rd entry of edges
function RS_CompareEdge( E1,E2 )
	if E1`r le E2`r then
		return 1;
	else
		return -1;
	end if;
end function;

// Weights for edges in spanning tree
function DE_Weight( P1,P2,P3 : Lambda := RPI/2 )
	z := (2*P3 - P2 - P1) / (P2 - P1);
	return Abs(Im(Argsinh(Argtanh(z)/Lambda)));
end function;
procedure DE_Edge_Weight( ~Edge, Points, Len : Lambda := RPI/2 )
	CCV := [];
	for k in [1..Len] do
		if k notin Edge`EP then
			uk := (2*Points[k] - Points[Edge`EP[2]] - Points[Edge`EP[1]]) / (Points[Edge`EP[2]] - Points[Edge`EP[1]]);
			Append(~CCV,uk);
			Edge`r := Min(Edge`r,Abs(Im(Argsinh(Argtanh(uk)/Lambda))));
		end if;
	end for;
	Edge`Data := CCV;
end procedure;
procedure GJ_Edge_Weight( ~Edge, Points, Len )
	V_r := [];
	for k in [1..Len] do
		if k notin Edge`EP then
			uk := (2*Points[k] - Points[Edge`EP[2]] - Points[Edge`EP[1]]) / (Points[Edge`EP[2]] - Points[Edge`EP[1]]);
			rk := (Abs(uk+1) + Abs(uk-1))/2;
			Append(~V_r,rk);
		end if;
	end for;
	Edge`Data := V_r;
	Edge`r := Min(Append(V_r,5.0));
end procedure;
// Weights for edges in Abel-Jacobi map
procedure DE_AJM_Weight( ~Edge, Points, Len : Lambda := RPI/2 )
	CCV := []; Append(~Edge,5.);
	for k in [1..Len] do
		if k ne Edge[3] then
			uk := (2*Points[k] - Edge[1][1] - Points[Edge[3]]) / (Edge[1][1] - Points[Edge[3]]);
			Append(~CCV,uk);
			Edge[4] := Min(Edge[4],Abs(Im(Argsinh(Argtanh(uk)/Lambda))));
		end if;
	end for;
	Append(~Edge,CCV);
end procedure;
procedure GJ_AJM_Weight( ~Edge, Points, Len )
	V_r := []; Append(~Edge,5.);
	for k in [1..Len] do
		if k ne Edge[3] then
			uk := (2*Points[k] - Edge[1][1] - Points[Edge[3]]) / (Edge[1][1] - Points[Edge[3]]);
			rk := (Abs(uk+1) + Abs(uk-1))/2;
			Append(~V_r,rk);
			Edge[4] := Min(Edge[4],rk);
		end if;
	end for;
	Append(~Edge,V_r);
end procedure;

procedure GJ_Params_Tree(STree,Points,m)
// Compute parameters for Gauss-Chebychev integration
	
	// Min_r too bad?
	assert STree`IntPars[1]-1 ge 3*10^-3;

	eps := 1/500;
	IntPars := [ STree`Edges[k]`r : k in [1..STree`Length ] ];
	rm := STree`IntPars[1];
	Groups := IntGroups(IntPars,rm,m);
	vprint SE,2 : "Groups:",Groups;
	if rm lt 1+(1/250) then
		Lr := [(1/2)*(rm+1)]; 
	else
		Lr := [ rm - eps];
	end if;
	Lr cat:= [ Min(g) - eps : g in Remove(Groups,1) | #g gt 1 ];	

	// Find best integration scheme for each edge and compute bound M
	NSchemes := #Lr;
	LrM := [ <1,Lr[l]> : l in [1..NSchemes] ];
	for k in [1..STree`Length] do
		M := 1;
		l := Max([ l : l in [1..NSchemes] | STree`Edges[k]`r gt Lr[l] ]);
		STree`Edges[k]`IntSch := l;
		M := Exp( -(1/m) * Log( &*[ STree`Edges[k]`Data[j] - Lr[l] : j in [1..STree`Length-1]]));
		if M ge 1 then
			M := M^(m-1);
		end if;
		M := Ceiling(Lr[l]^(STree`Length-1) * M);
		LrM[l][1] := Max(M,LrM[l][1]);
	end for;

	vprint SE,2 : "Parameters(tree):",LrM;
	// Save parameters
	STree`Params := LrM;
end procedure;

function GJ_Params_AJM(ComplexEdges,SEC)
// Compute parameters for Gauss-Chebychev integration
	
	// Edges = <[P],v_P,k,r_E>
	m := SEC`Degree[1]; n := SEC`Degree[2];
	NEdges := #ComplexEdges; NewEdges := [];
	
	rm := 5.;
	// NewEdges = <[P],v_P,k,r_E,V_r>
	for j in [1..NEdges] do
		Edge := ComplexEdges[j];
		GJ_AJM_Weight(~Edge,SEC`LowPrecBranchPoints,n);
		rm := Min(rm,Edge[4]);
		Append(~NewEdges,Edge);
	end for;

	eps := 1/500;
	IntPars := [ NewEdges[k][4] : k in [1..NEdges] ];
	Groups := IntGroups(IntPars,rm,m);
	vprint SE,2 : "Groups:",Groups;
	if rm lt 1+(1/250) then
		Lr := [(1/2)*(rm+1)]; 
	else
		Lr := [ rm - eps];
	end if;
	Lr cat:= [ Min(g) - eps : g in Remove(Groups,1) | #g gt 1 ];

	// Find best integration scheme for each edge and compute bound M
	NSchemes := #Lr;
	LrM := [ <1,Lr[l]> : l in [1..NSchemes] ];
	for k in [1..#ComplexEdges] do
		M := 1;
		l := Max([ l : l in [1..NSchemes] | NewEdges[k][4] gt Lr[l] ]);
		NewEdges[k][4] := l;
		M := Exp( -(1/m) * Log( &*[ NewEdges[k][5][j] - Lr[l] : j in [1..SEC`Degree[2]-1]]));
		if M ge 1 then
			M := M^(m-1);
		end if;
		M := Ceiling(Lr[l]^(SEC`Degree[2]-2) * M);
		LrM[l][1] := Max(M,LrM[l][1]);
	end for;

	vprint SE,2 : "Parameters(AJM):",LrM;

	// Return parameters
	return LrM, NewEdges;
end function;

procedure DE_Params_Tree(STree,Points,m)
// Computes double-exponential integration integration parameters for a spanning tree
	
	// Make list of r's
	MaxMinDiff := STree`IntPars[2]-STree`IntPars[1];
	NSchemes := Max(Min(Floor(STree`Length/2),Floor(10*MaxMinDiff)),1);

	vprint SE,2 : "Number of schemes:",NSchemes;
	r_min := (19/20) * STree`IntPars[1];
	r_max := (19/20) * STree`IntPars[2];

	if NSchemes eq 1 and Abs(r_min-r_max) lt 10^-5 then
		Lr := [ r_min ];
	else
		Lr := [ (1-t/NSchemes)*r_min + t/NSchemes*r_max : t in [0..NSchemes] ];
	end if;

	assert &and[STree`IntPars[2] lt RPI, STree`IntPars[1] gt 0];

	// Find best integration scheme for each edge and compute bounds M1,M2
	NSchemes := #Lr; 
	LrM := [ <1,1,Lr[l]> : l in [1..NSchemes] ];
	for k in [1..STree`Length] do
		l := NSchemes;
		while (19/20)*STree`Edges[k]`r lt Lr[l] do
			l -:= 1;
		end while;
		STree`Edges[k]`IntSch := l;
		M1 := Ceiling(Bound_M1(STree`Edges[k]`Data,STree`Length-1,m));
		M2 := Ceiling(Bound_M2(STree`Edges[k]`Data,STree`Length-1,m,STree`Length-1,Lr[l]));
		assert M2 ge M1;
		LrM[l][1] := Max(M1,LrM[l][1]);
		LrM[l][2] := Max(M2,LrM[l][2]);
	end for;
	vprint SE,2 : "Parameters(tree):",LrM;
	// Save parameters
	STree`Params := LrM;
end procedure;


function DE_Params_AJM(ComplexEdges,SEC)
// Computes double-exponential integration integration parameters for a spanning tree

	// Edges = <[P],v_P,k,r_E>
	m := SEC`Degree[1]; 
	n := SEC`Degree[2];
	NEdges := #ComplexEdges;
	NewEdges := [];
	
	Min_r := 5.;
	Max_r := 0.;	

	for j in [1..NEdges] do
		Edge := ComplexEdges[j];
		DE_AJM_Weight(~Edge,SEC`LowPrecBranchPoints,n);
		Min_r := Min(Min_r,Edge[4]);
		Max_r := Max(Max_r,Edge[4]);
		Append(~NewEdges,Edge);
	end for;
	// NewEdges = <[P],v_P,k,r_E,CCV>

	// Make list of r's
	Min_r *:= (9/10);
	Max_r *:= (9/10);
	
	NSchemes := Max(Min(Floor(NEdges/2),Floor(10*(Max_r-Min_r))),1);

	vprint SE,2 : "Number of schemes(AJM):",NSchemes;
	if NSchemes eq 1 and Abs(Min_r-Max_r) lt 10^-5 then
		Lr := [ Min_r ];
	else
		Lr := [ (1-t/NSchemes)*Min_r + t/NSchemes*Max_r : t in [0..NSchemes] ];
	end if;
	assert &and[Max_r lt RPI, Min_r gt 0];

	// Find best integration scheme for each edge
	NSchemes := #Lr; 
	LrM := [ <1,1,Lr[l]> : l in [1..NSchemes] ];
	for k in [1..NEdges] do
		l := NSchemes;
		while (19/20)*NewEdges[k][4] lt Lr[l] do
			l -:= 1;
		end while;
		NewEdges[k][4] := l;
		M1 := Ceiling(Bound_M1(NewEdges[k][5],n-1,m:AJM));
		M2 := Ceiling(Bound_M2(NewEdges[k][5],n-1,m,n-2,Lr[l]:AJM));
		LrM[l][1] := Max(M1,LrM[l][1]);
		LrM[l][2] := Max(M2,LrM[l][2]);
	end for;

	vprint SE,2 : "Parameters(AJM):",LrM;

	// Return parameters
	return LrM, NewEdges;
end function;


// Edges
declare type SEEdge;

declare attributes SEEdge :
	EP, // End points = [E1,E2]
	Data, // Data = [ u_1, ... , u_{n-2},(b-a)/2,(b+a)/(b-a),C_ab ]
	up,
	r,
	IntSch;

intrinsic Print( E::SEEdge )
{ Printing }
	print "Endpoints:",E`EP;
	/*print "Data:",E`Data;
	print "up:",E`up;
	print "r:",E`r;
	print "IntSch:",E`IntSch;*/
end intrinsic;


function SE_Edge(k,l:r:=0)
	E := New(SEEdge);
	E`EP := [k,l];
	if r ne 0 then
		E`r := r;
	end if;
	return E;
end function;

procedure EdgeData( E,Points,Zetas,m,n )
	a := Points[E`EP[1]];
	b := Points[E`EP[2]];
	if E`EP[1] lt E`EP[2] then
		Pts := Remove(Remove(Points,E`EP[1]),E`EP[2]-1);
	else
		Pts := Remove(Remove(Points,E`EP[2]),E`EP[1]-1);
	end if;
	CCV, up := SortByRealPart([ (2*x-b-a)/(b-a) : x in Pts ]);
	bma := b-a;
	Append(~CCV,bma/2);
	Append(~CCV,(b+a)/bma);
	C_ab := Exp( (n/m) * Log(b-a)) * Zetas[(up+1) mod 2 + 1]; 
	Append(~CCV,C_ab);
	E`up := up;
	E`Data := CCV;
end procedure;

procedure FlipEdge( E )
	E`EP := Reverse(E`EP);
end procedure;


procedure TreeData( ~SpTree, Points, Zetas, m )
	//SpTree`Data := [ EdgeData(E,Points) : E in SpTree`Edges ];
	for E in SpTree`Edges do
		EdgeData(E,Points,Zetas,m,SpTree`Length+1);
	end for;
end procedure;

// Spanning tree type 
declare type SpTree;

declare attributes SpTree :
	Length,
	Edges,
	IntPars,
	ExtEdges,
	Params;

procedure SpanningTree(SEC)
// Computes a spanning tree between Points

	assert SEC`IntegrationType in ["DE","GJ"];
	Points := SEC`LowPrecBranchPoints;	
	Len := SEC`Degree[2];
	Edges := [];

	T := New(SpTree);
	T`Length := Len-1;
	T`Edges := [];
	T`ExtEdges := [];

	Min_r := 5.;
	Max_r := 0.;	

	// Use minimal euclidean distance tree?
	if true then
	//if 2*Len gt SEC`Prec then
		EuclideanDistance := true;
	else
		EuclideanDistance := false;
	end if;

	for k in [1..Len] do
		for l in [k+1..Len] do
			if EuclideanDistance then
				//Edge := <k,l,-Abs(Points[k]-Points[l])>;
				Edge := SE_Edge(k,l:r:=-Abs(Points[k]-Points[l]));
			else
				//Edge := <k,l,20.>;
				Edge := SE_Edge(k,l:r:=5.);
				if SEC`IntegrationType eq "DE" then
					DE_Edge_Weight(~Edge,Points,Len);
				elif SEC`IntegrationType eq "GJ" then
					GJ_Edge_Weight(~Edge,Points,Len);
				else
					error Error("Unsupported integration type.");
				end if;
			end if;
			Append(~Edges,Edge);
		end for;
	end for;
	
	// Sort by third entry
	Sort(~Edges,RS_CompareEdge);

	Taken := [ 0 : j in [1..Len] ];
	k := 0;
	while k lt T`Length do
		l := 1;
		if k ne 0 then
			while Taken[Edges[l]`EP[1]] eq Taken[Edges[l]`EP[2]] do
				l +:= 1;
			end while;
		end if;
		if EuclideanDistance then
			NewEdge := SE_Edge(Edges[l]`EP[1],Edges[l]`EP[2]:r:=5.);
			if SEC`IntegrationType eq "DE" then
				DE_Edge_Weight(~NewEdge,Points,Len);
			elif SEC`IntegrationType eq "GJ" then
				GJ_Edge_Weight(~NewEdge,Points,Len);
			else
				error Error("Unsupported integration type.");
			end if;
		else
			NewEdge := Edges[l];
		end if;
		if Taken[NewEdge`EP[2]] eq 1 then
			// Flip edge
			FlipEdge(NewEdge);
		end if;
		Append(~T`Edges,NewEdge);
		k +:= 1;
		Taken[NewEdge`EP[1]] := 1;
		Taken[NewEdge`EP[2]] := 1;
		
		if Min_r gt NewEdge`r then
			Min_r := NewEdge`r;
			T`ExtEdges[1] := k;
		end if;	
		if Max_r lt NewEdge`r then
			Max_r := NewEdge`r;
			T`ExtEdges[2] := k;
		end if;
	end while;
	T`IntPars := [Min_r,Max_r];

	if SEC`IntegrationType eq "GJ" then
		if T`IntPars[1]-1 lt 3*10^-3 then
			SEC`IntegrationType := "DE";
			vprint SE,1 : "Changed integration type to",SEC`IntegrationType,"due to bad integration parameter:",T`IntPars[1];
			SpanningTree(SEC);
		else
			SEC`SpanningTree := T;
			GJ_Params_Tree(SEC`SpanningTree,Points,SEC`Degree[1]);
		end if;
	else
		SEC`SpanningTree := T;
		DE_Params_Tree(SEC`SpanningTree,Points,SEC`Degree[1]);		
	end if;

end procedure;

intrinsic Print( STree::SpTree : Edges := false )
{ Printing }
	print "Spanning tree between",STree`Length+1,"points";
	print "Integration parameter:",STree`IntPars;
	print "Worst/Best edge:",STree`ExtEdges;
	//print "Params:",STree`Params;
	if Edges then
		print "Edges:",STree`Edges;
	end if;
end intrinsic;



