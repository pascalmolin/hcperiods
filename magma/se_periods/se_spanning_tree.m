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
function GC_Weight( P1,P2,P3 )
	return (Abs(P3-P1)+Abs(P3-P2))/Abs(P2-P1);
end function;
procedure GC_Edge_Weight( ~Edge, Points, Len )
	for k in [1..Len] do
		if k notin [Edge[1],Edge[2]] then
			Edge[3] := Min(Edge[3],GC_Weight(Points[Edge[1]],Points[Edge[2]],Points[k]));
		end if;
	end for;
end procedure;
// Weights for edges in Abel-Jacobi map
function DE_AJM_Weight( E, Points, Len )
	r := 4.;
	for k in [1..Len] do
		if k ne E[1] then
			r := Min(r,DE_Weight(E[2],Points[E[1]],Points[k]));
		end if;
	end for;
	return r;
end function;;
// Weights for edges in Abel-Jacobi map
function GC_AJM_Weight( E, Points, Len )
	r := 4.;
	for k in [1..Len] do
		if k ne E[1] then
			r := Min(r,GC_Weight(E[2],Points[E[1]],Points[k]));
		end if;
	end for;
	return r;
end function;;

procedure GC_Params_Tree(STree,Points,Prec)
// Compute parameters for Gauss-Chebychev integration
	
	// Compute data for worst edge
	E := STree`Edges[STree`WorstEdge];
	CCV, up := MakeCCVector( E,Points );

	// Compute r_k and r_0
	V_r := []; r_0 := 4.;
	for k in [1..STree`Length-1] do
		r_k := (1/2) * ( Abs(CCV[k]+1) + Abs(CCV[k]-1) );
		assert r_k gt 1;
		r_0 := Min(r_0,r_k);
		Append(~V_r,r_k);
	end for;
	vprint SE,2 : "r_0:",r_0;
	vprint SE,3 : "V_r:",V_r;

   	// Make it precise
	r := (1/2) * (r_0 + 1);
	rf := Floor(r_0);
	assert rf in [1..4];
	vprint SE,3 : "rf:",rf;
	r := rf/(rf+1)*r_0 + 1/(rf+1);
	

	//r := (9/10) * r_0 + (1/10);
	g := Floor((STree`Length)/2);
	M_r := r^(g-1)/Sqrt(&*[V_r[k] - r : k in [1..STree`Length-1]]);
	M_r := Max(1,Ceiling(M_r));

	vprint SE,2 : "r:",r;
	vprint SE,2 : "M(r):",M_r;

	// Number of integration points
	N := Ceiling((Log(10)*Prec+Log(2*Pi20*M_r)+1)/(2*Argcosh(r)));

	STree`Params := < M_r, N, r>;
	vprint SE,2 : "N:",N;
end procedure;


procedure DE_Params_Tree(STree,Points,m,Prec)
// Computes double-exponential integration integration parameters for a spanning tree
	
	// Check r
	r := STree`IntPar;
	if r ge Pi20/2 then
		error Error("Not supposed to happen.");
	end if;
	assert &and[r lt Pi20, r gt 0];

	// Make it precise
	r *:= (19/20);	
	//r *:= 9/10;
	vprint SE,2 : "r:",r;

	// Compute data for worst edge
	E := STree`Edges[STree`WorstEdge];
	CCV, up := MakeCCVector(E,Points);

	// Compute bounds M1,M2 for worst edge
	M1 := Ceiling(Bound_M1(CCV,STree`Length-1,m));
	M2 := Ceiling(Bound_M2(CCV,STree`Length-1,m,r));
	vprint SE,2 : "M1:",M1;
	vprint SE,2 : "M2:",M2;

	// Save parameters
	STree`Params := <M1, M2, r>;
end procedure;

function DE_Params_AJM(Edges,SEC)
// Computes double-exponential integration integration parameters for a spanning tree

	// Edges = <[P],v_P,k>
	m := SEC`Degree[1]; n := SEC`Degree[2];
	r := 5.;
	// Find worst edge and compute r
	for E in Edges do
		r_E := DE_AJM_Weight(<E[3],E[1][1]>,SEC`LowPrecBranchPoints,n);
		if r_E lt r then
			r := r_E;
			MinEdge := E;
		end if;
	end for;
	if r ge Pi20/2 then
		error Error("Not supposed to happen.");
	end if;
	assert &and[r lt Pi20/2, r gt 0];
	
	// Make it precise
	r *:= (19/20);	
	//r *:= (9/10);

	// Compute bounds M1,M2
	LowPrecData := MakeCCVector(<MinEdge[3],MinEdge[1][1]>,SEC`LowPrecBranchPoints);
	M1 := Ceiling(Bound_M1(LowPrecData,n-1,m));
	M2 := Ceiling(Bound_M2(LowPrecData,n-1,m,r));

	// Return parameters
	return <M1, M2, r>;
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
	IntPar,
	WorstEdge,
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

	if SEC`Degree[2] ge 15 then
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

	T`IntPar := 4.;
	Taken := [ 0 : j in [1..Len] ];
	k := 0;
	while k lt T`Length do
		l := 1;
		if k ne 0 then
			while Taken[Edges[l][1]] eq Taken[Edges[l][2]] do
				l +:= 1;
			end while;
		end if;

		if Taken[Edges[l][2]] eq 1 then
			// Flip edge
			Append(~T`Edges,<Edges[l][2],Edges[l][1]>);
		else
			Append(~T`Edges,<Edges[l][1],Edges[l][2]>);	
		end if;
		Taken[Edges[l][1]] := 1;
		Taken[Edges[l][2]] := 1;
		if ED then
			Edges[l][3] := 4.;
			if SEC`IntegrationType eq "DE" then
				DE_Edge_Weight(~Edges[l],Points,Len);
			elif SEC`IntegrationType eq "GC" then
				GC_Edge_Weight(~Edges[l],Points,Len);
			else
				error Error("Unsupported integration type.");
			end if;
		end if;
		k +:= 1;
		if T`IntPar gt Edges[l][3] then
			T`IntPar := Edges[l][3];
			T`WorstEdge := k;
		end if;	
	end while;

	if SEC`IntegrationType eq "GC" and T`IntPar-1 lt 10^-3 then
		SEC`IntegrationType := "DE";
		vprint SE,1 : "Changed integration type to",SEC`IntegrationType,"due to bad integration parameter:",T`IntPar;
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
	print "Integration parameter:",STree`IntPar;
	print "Worst edge:",STree`Edges[STree`WorstEdge];
	print "Params:",STree`Params;
	print "Edges:",STree`Edges;
end intrinsic;



