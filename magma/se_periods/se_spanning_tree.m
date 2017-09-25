/*******************************************************************************

 Copyright (C) 2017 
 Code adapted from Pascal Molin,
 written by Christian Neurohr

*******************************************************************************/

// Import functions
import "se_help_funcs.m": MakeCCVector;
import "se_de_int.m": Bound_M1, Bound_M2;
import "se_gc_int.m": DistanceEllipse;


// Complex field
C<I> := ComplexField(20);
RPI := Real(Pi(C));

// Compare 3rd entry of edges
function RS_CompareEdge( E1,E2 )
	if E1[3] le E2[3] then
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
		if k notin [Edge[1],Edge[2]] then
			uk := (2*Points[k] - Points[Edge[2]] - Points[Edge[1]]) / (Points[Edge[2]] - Points[Edge[1]]);
			Append(~CCV,uk);
			Edge[3] := Min(Edge[3],Abs(Im(Argsinh(Argtanh(uk)/Lambda))));
		end if;
	end for;
	Append(~Edge,CCV);
end procedure;
procedure GJ_Edge_Weight( ~Edge, Points, Len )
	V_r := [];
	for k in [1..Len] do
		if k notin [Edge[1],Edge[2]] then
			uk := (2*Points[k] - Points[Edge[2]] - Points[Edge[1]]) / (Points[Edge[2]] - Points[Edge[1]]);
			rk := (Abs(uk+1) + Abs(uk-1))/2;
			Append(~V_r,rk);
		end if;
	end for;
	Append(~Edge,V_r);
	Edge[3] := Min(V_r);
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


procedure GJ_Params_Tree(STree,Points,m,Prec)
// Compute parameters for Gauss-Chebychev integration
	
	// Min_r too bad?
	assert STree`IntPars[1]-1 ge 3*10^-3;

	OT := 1/1000;
	Lr := []; r := STree`IntPars[1];
	if r lt 1+(1/500) then
		r := (1/2)*(r+1); 
	else
		r -:= OT;
	end if;
	if m eq 2 then
		k := 1/8;
	else
		k := 10;
	end if;
	while r lt STree`IntPars[2] do
		Append(~Lr,r);
		r +:= (k/10);
		k *:= 2;
	end while;

	// Find best integration scheme for each edge and compute bound M
	gm1 := Floor((STree`Length)/2);
	NSchemes := #Lr;
	LrM := [ <1,Lr[l]> : l in [1..NSchemes] ];
	for k in [1..STree`Length] do
		M := 1;
		l := NSchemes;
		while STree`Edges[k][3]-OT lt Lr[l] do
			l -:= 1;
		end while;
		STree`Edges[k][3] := l;
		M := Ceiling(Lr[l]^gm1 * Exp( (1/m) * Log( &*[ STree`Edges[k][4][j] - Lr[l] : j in [1..STree`Length-1]])));
		LrM[l][1] := Max(M,LrM[l][1]);
	end for;

	// Number of integration points N
	LrMN := [];
	for k in [1..NSchemes] do
		P := LrM[k];
		achr := Argcosh(P[2]);
		//Append(~LrMN,Append(P,Ceiling((Log(10)*Prec+Log(2*RPI*P[1])+1)/(2*achr))));
		Append(~LrMN,Append(P,Ceiling((Log(32*P[1]/15)+Log(10)*Prec-Log(1-Exp(achr)^(-2)))/(2*achr))));
	end for;
	vprint SE,2 : "Parameters(tree):",LrMN;
	// Save parameters
	STree`Params := LrMN;
end procedure;


procedure DE_Params_Tree(STree,Points,m,Prec)
// Computes double-exponential integration integration parameters for a spanning tree
	
	// Make list of r's
	MaxMinDiff := STree`IntPars[2]-STree`IntPars[1];
	NSchemes := Max(Min(STree`Length,Floor(20*MaxMinDiff)),1);

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
		while (19/20)*STree`Edges[k][3] lt Lr[l] do
			l -:= 1;
		end while;
		STree`Edges[k][3] := l;
		M1 := Ceiling(Bound_M1(STree`Edges[k][4],STree`Length-1,m));
		M2 := Ceiling(Bound_M2(STree`Edges[k][4],STree`Length-1,m,STree`Length-1,Lr[l]));
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
	m := SEC`Degree[1]; n := SEC`Degree[2];
	NEdges := #ComplexEdges; NewEdges := [];
	
	Min_r := 5.;
	Max_r := 0.;	

	for j in [1..NEdges] do
		Edge := ComplexEdges[j];
		DE_AJM_Weight(~Edge,SEC`LowPrecBranchPoints,n);
		if Edge[4] lt Min_r then
			Min_r := Edge[4];
			MinEdge := Edge;
		end if;
		Max_r := Max(Max_r,Edge[4]);
		Append(~NewEdges,Edge);
	end for;
	// NewEdges = <[P],v_P,k,r_E,CCV>

	// Make list of r's
	Min_r *:= (19/20);
	Max_r *:= (19/20);
	
	NSchemes := Max(Min(NEdges,Floor(20*(Max_r-Min_r))),1);

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
		assert M2 ge M1;
		LrM[l][1] := Max(M1,LrM[l][1]);
		LrM[l][2] := Max(M2,LrM[l][2]);
	end for;

	vprint SE,2 : "Parameters(AJM):",LrM;

	// Return parameters
	return LrM, NewEdges;
end function;


function EdgeData( E,Points )
	// Data = [ u_1, ... , u_{n-2},(b-a)/2,(b+a)/(b-a),up ]
	Data, up :=  MakeCCVector( E,Points );
	return Append(Data,up);
end function;

procedure TreeData( ~SpTree, Points )
	SpTree`Data := [ EdgeData(E,Points) : E in SpTree`Edges ];
end procedure;

// Spanning tree
declare type SpTree;

declare attributes SpTree :
	Length,
	Edges,
	IntPars,
	ExtEdges,
	Params,
	Data;


procedure SpanningTree(SEC)
// Computes a spanning tree between Points

	assert SEC`IntegrationType in ["DE","GC","GJ"];
	Points := SEC`LowPrecBranchPoints;	
	Len := SEC`Degree[2];
	Edges := [];

	T := New(SpTree);
	T`Length := Len-1;
	T`Edges := [];
	T`ExtEdges := [];

	Min_r := 20.;
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
				Edge := <k,l,-Abs(Points[k]-Points[l])>;
			else
				Edge := <k,l,20.>;
				if SEC`IntegrationType eq "DE" then
					DE_Edge_Weight(~Edge,Points,Len);
				elif SEC`IntegrationType in ["GC","GJ"] then
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
			while Taken[Edges[l][1]] eq Taken[Edges[l][2]] do
				l +:= 1;
			end while;
		end if;
		if EuclideanDistance then
			NewEdge := <Edges[l][1],Edges[l][2],20.>;
			if SEC`IntegrationType eq "DE" then
				DE_Edge_Weight(~NewEdge,Points,Len);
			elif SEC`IntegrationType in ["GC","GJ"] then
				GJ_Edge_Weight(~NewEdge,Points,Len);
			else
				error Error("Unsupported integration type.");
			end if;
		else
			NewEdge := Edges[l];
		end if;
		if Taken[NewEdge[2]] eq 1 then
			// Flip edge
			Append(~T`Edges,<NewEdge[2],NewEdge[1],NewEdge[3],NewEdge[4]>);
		else
			Append(~T`Edges,NewEdge);	
		end if;
		k +:= 1;
		Taken[NewEdge[1]] := 1;
		Taken[NewEdge[2]] := 1;
		
		if Min_r gt NewEdge[3] then
			Min_r := NewEdge[3];
			T`ExtEdges[1] := k;
		end if;	
		if Max_r lt NewEdge[3] then
			Max_r := NewEdge[3];
			T`ExtEdges[2] := k;
		end if;
	end while;
	T`IntPars := [Min_r,Max_r];

	if SEC`IntegrationType in ["GC","GJ"] then
		if T`IntPars[1]-1 lt 3*10^-3 then
			SEC`IntegrationType := "DE";
			vprint SE,1 : "Changed integration type to",SEC`IntegrationType,"due to bad integration parameter:",T`IntPars[1];
			SpanningTree(SEC);
		else
			SEC`SpanningTree := T;
			GJ_Params_Tree(SEC`SpanningTree,Points,SEC`Degree[1],SEC`Prec);
		end if;
	else
		SEC`SpanningTree := T;
		DE_Params_Tree(SEC`SpanningTree,Points,SEC`Degree[1],SEC`Prec);		
	end if;

end procedure;

intrinsic Print( STree::SpTree : Edges := false )
{ Printing }
	print "Spanning tree between",STree`Length+1,"points";
	print "Integration parameter:",STree`IntPars;
	print "Worst/Best edge:",STree`ExtEdges;
	print "Params:",STree`Params;
	if Edges then
		print "Edges:",STree`Edges;
	end if;
end intrinsic;



