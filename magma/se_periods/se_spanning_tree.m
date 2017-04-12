/******************************************************************************

 Copyright (C) 2017 Christian Neurohr

 ******************************************************************************/

// Import functions
import "se_help_funcs.m": MakeCCVector;


// Complex field
C20<i> := ComplexField(20);
Pi20 := Pi(C20);

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
		print "E:",E;
		if k ne E[1] then
			print "E[1]:",E[1];
			print "t:",Type(E[1]);
			print "Points:",Points;
			r := Min(r,DE_Weight(E[2],Points[E[1]],Points[k]));
		end if;
	end for;
	return r;
end function;;


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
	Len,
	Edges,
	IntPar,
	WorstEdge,
	Data;


function SpanningTree(S)
// Computes a spanning tree between Points

	assert S`IntegrationType in ["DE","GC"];
	Points := S`LowPrecBranchPoints;	
	Len := S`Degree[2];
	Edges := [];

	T := New(SpTree);
	T`Len := Len-1;
	T`Edges := [];

	if S`Degree[2] ge 15 then
		ED := true;
	else
		ED := false;
	end if;

	for k in [1..Len] do
		for l in [k+1..Len] do
			if ED then
				Edge := <k,l,-Abs(Points[k]-Points[l])>;
			else
				Edge := <k,l,4.>;
				if S`IntegrationType eq "DE" then
					DE_Edge_Weight(~Edge,Points,Len);
				elif S`IntegrationType eq "GC" then
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
	while k lt T`Len do
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
			if S`IntegrationType eq "DE" then
				DE_Edge_Weight(~Edges[l],Points,Len);
			elif S`IntegrationType eq "GC" then
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
	return T;
end function;

intrinsic Print( STree::SpTree )
{ Printing }
	print "Spanning tree between",STree`Len+1,"points";
	print "Integration parameter:",STree`IntPar;
	print "Worst edge:",STree`Edges[STree`WorstEdge];
	print "Edges:",STree`Edges;
end intrinsic;



