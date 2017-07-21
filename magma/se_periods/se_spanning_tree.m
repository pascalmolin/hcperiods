/*******************************************************************************

 Copyright (C) 2017 
 Code adapted from Pascal Molin,
 written by Christian Neurohr

*******************************************************************************/

// Import functions
import "se_help_funcs.m": MakeCCVector;
import "se_de_int.m": Bound_M1, Bound_M2;


// Complex field
C20<I> := ComplexField(20);
Pi20 := Real(Pi(C20));

// Compare 3rd entry of edges
function RS_CompareEdge( E1,E2 )
	if E1[3] le E2[3] then
		return 1;
	else
		return -1;
	end if;
end function;

// Weights for edges in spanning tree
function DE_Weight( P1,P2,P3 : Lambda := Pi20/2 )
	z := (2*P3 - P2 - P1) / (P2 - P1);
	return Abs(Im(Argsinh(Argtanh(z)/Lambda)));
end function;
procedure DE_Edge_Weight( ~Edge, Points, Len )
	for k in [1..Len] do
		if k notin [Edge[1],Edge[2]] then
			Edge[3] := Min(Edge[3],DE_Weight(Points[Edge[1]],Points[Edge[2]],Points[k]));
		end if;
	end for;
end procedure;
procedure GC_Edge_Weight( ~Edge, Points, Len )
	CCV, up := MakeCCVector( Edge,Points );
	r_0 := 5.;
	for k in [1..Len-2] do
		r_k := (1/2) * ( Abs(CCV[k]+1) + Abs(CCV[k]-1) );
		r_0 := Min(r_0,r_k);
	end for;
	Edge[3] := r_0;
end procedure;
// Weights for edges in Abel-Jacobi map
procedure DE_AJM_Weight( ~Edge, Points, Len )
	r := 5.;
	for k in [1..Len] do
		if k ne Edge[3] then
			r := Min(r,DE_Weight(Edge[1][1],Points[Edge[3]],Points[k]));
		end if;
	end for;
	Append(~Edge,r);
end procedure;


procedure GC_Params_Tree(STree,Points,Prec)
// Compute parameters for Gauss-Chebychev integration
	
	// Make list of r's
	MaxMinDiff := STree`IntPars[2]-STree`IntPars[1];
	NSchemes := Max(Min(STree`Length,Floor(3*MaxMinDiff)),1);

	// Min_r too bad?
	assert STree`IntPars[1]-1 gt 10^-3;

   	vprint SE,2 : "Number of schemes:",NSchemes;
	if NSchemes eq 1 and Abs(STree`IntPars[2]-STree`IntPars[2]) lt 10^-10 then
		Lr := [ STree`IntPars[1] - (1/1000) ];
	else
		Avg_r := (1/2)*(STree`IntPars[2] + STree`IntPars[1]);
		Lr := [ (1-t/NSchemes)*STree`IntPars[1] + t/NSchemes*Avg_r - (1/1000) : t in [0..NSchemes] ];
	end if;
	vprint SE,2 : "Lr:",Lr;

	// Compute data for worst edge
	E := STree`Edges[STree`ExtEdges[1]];
	CCV, up := MakeCCVector( E,Points );
	
	// Compute bound M(r)
	V_r := [ (1/2) * ( Abs(CCV[k]+1) + Abs(CCV[k]-1)) : k in [1..STree`Length-1] ];
	g := Floor((STree`Length)/2);
	M_r := Lr[1]^(g-1)/Sqrt(&*[V_r[k] - Lr[1] : k in [1..STree`Length-1]]);
	M_r := Max(1,Ceiling(M_r));

	vprint SE,2 : "r:",Lr[1];
	vprint SE,2 : "M(r):",M_r;

	// Number of integration points N
	NPoints := [];
	for r in Lr do
		achr := Argcosh(r);
		//N := Ceiling((Log(10)*Prec+Log(2*Pi20*M_r)+1)/(2*achr));
		Append(~NPoints,Ceiling((Log(64*M_r/15)+Log(10)*Prec-Log(1-Exp(achr)^(-2)))/(2*achr)));
	end for;
	STree`Params := <M_r, NPoints, Lr>;
	vprint SE,2 : "N:",NPoints;
end procedure;


procedure DE_Params_Tree(STree,Points,m,Prec)
// Computes double-exponential integration integration parameters for a spanning tree
	
	// Make list of r's
	MaxMinDiff := STree`IntPars[2]-STree`IntPars[1];
	NSchemes := Max(Min(STree`Length,Floor(20*MaxMinDiff)),1);

	vprint SE,2 : "Number of schemes:",NSchemes;
	if NSchemes eq 1 and Abs(STree`IntPars[2]-STree`IntPars[2]) lt 10^-10 then
		Lr := [ (19/20) * STree`IntPars[1] ];
	else
		Lr := [ (19/20) * ((1-t/NSchemes)*STree`IntPars[1] + t/NSchemes*STree`IntPars[2]) : t in [0..NSchemes] ];
	end if;
	vprint SE,2 : "Lr:",Lr;

	assert &and[STree`IntPars[2] lt Pi20, STree`IntPars[1] gt 0];

	// Compute data for worst edge
	E := STree`Edges[STree`ExtEdges[1]];
	CCV, up := MakeCCVector(E,Points);

	// Compute bounds M1,M2 for worst edge
	M1 := Ceiling(Bound_M1(CCV,STree`Length-1,m));
	M2 := Ceiling(Bound_M2(CCV,STree`Length-1,m,Lr[1]));
	vprint SE,2 : "M1:",M1;
	vprint SE,2 : "M2:",M2;

	// Save parameters
	STree`Params := <M1, M2, Lr>;
end procedure;


function DE_Params_AJM(ComplexEdges,SEC)
// Computes double-exponential integration integration parameters for a spanning tree

	// Edges = <[P],v_P,k,r_E>
	m := SEC`Degree[1]; n := SEC`Degree[2]; NEdges := #ComplexEdges; NewEdges := [];
	
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

	// Make list of r's
	MaxMinDiff := Max_r-Min_r;
	NSchemes := Max(Min(NEdges,Floor(20*MaxMinDiff)),1);

	vprint SE,2 : "Number of schemes:",NSchemes;
	if NSchemes eq 1 and Abs(Min_r-Max_r) lt 10^-10 then
		Lr := [ (19/20) * Min_r ];
	else
		Lr := [ (19/20) * ((1-t/NSchemes)*Min_r + t/NSchemes*Max_r) : t in [0..NSchemes] ];
	end if;

	assert &and[Max_r lt Pi20, Min_r gt 0];

	// Compute bounds M1,M2
	LowPrecData := MakeCCVector(<MinEdge[3],MinEdge[1][1]>,SEC`LowPrecBranchPoints);
	M1 := Ceiling(Bound_M1(LowPrecData,n-1,m));
	M2 := Ceiling(Bound_M2(LowPrecData,n-1,m,Min_r));

	// Return parameters
	return <M1, M2, Lr>, NewEdges;
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

	assert SEC`IntegrationType in ["DE","GC"];
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
	if 2*Len gt SEC`Prec then
		ED := true;
	else
		ED := false;
	end if;

	for k in [1..Len] do
		for l in [k+1..Len] do
			if ED then
				Edge := <k,l,-Abs(Points[k]-Points[l])>;
			else
				Edge := <k,l,5.>;
				if SEC`IntegrationType eq "DE" then
					DE_Edge_Weight(~Edge,Points,Len);
				elif SEC`IntegrationType eq "GC" then
					GC_Edge_Weight(~Edge,Points,Len);
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
		if ED then
			Edges[l][3] := 5.;
			if SEC`IntegrationType eq "DE" then
				DE_Edge_Weight(~Edges[l],Points,Len);
			elif SEC`IntegrationType eq "GC" then
				GC_Edge_Weight(~Edges[l],Points,Len);
			else
				error Error("Unsupported integration type.");
			end if;
		end if;
		if Taken[Edges[l][2]] eq 1 then
			// Flip edge
			Append(~T`Edges,<Edges[l][2],Edges[l][1],Edges[l][3]>);
		else
			Append(~T`Edges,Edges[l]);	
		end if;
		k +:= 1;
		Taken[Edges[l][1]] := 1;
		Taken[Edges[l][2]] := 1;
		
		if Min_r gt Edges[l][3] then
			Min_r := Edges[l][3];
			T`ExtEdges[1] := k;
		end if;	
		if Max_r lt Edges[l][3] then
			Max_r := Edges[l][3];
			T`ExtEdges[2] := k;
		end if;
	end while;
	T`IntPars := [Min_r,Max_r];

	if SEC`IntegrationType eq "GC" and T`IntPars[1]-1 lt 3*10^-3 then
		SEC`IntegrationType := "DE";
		vprint SE,1 : "Changed integration type to",SEC`IntegrationType,"due to bad integration parameter:",T`IntPars[1];
		SpanningTree(SEC);
	else
		SEC`SpanningTree := T;
	end if;

	if SEC`IntegrationType eq "GC" then
		GC_Params_Tree(SEC`SpanningTree,Points,SEC`Prec);
	else
		DE_Params_Tree(SEC`SpanningTree,Points,SEC`Degree[1],SEC`Prec);
	end if;
end procedure;

intrinsic Print( STree::SpTree )
{ Printing }
	print "Spanning tree between",STree`Length+1,"points";
	print "Integration parameter:",STree`IntPars;
	print "Worst/Best edge:",STree`ExtEdges;
	print "Params:",STree`Params;
	print "Edges:",STree`Edges;
end intrinsic;



