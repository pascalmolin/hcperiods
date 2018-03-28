/*******************************************************************************

 Copyright (C) 2017 
 Code adapted from Pascal Molin,
 written by Christian Neurohr

*******************************************************************************/

// Import functions
import "se_help_funcs.m": SortByRealPart, IntGroups;
import "se_de_int.m": Bound_M1, Bound_M2;
import "se_gj_int.m": DistanceEllipse;
import "se_anal_cont.m": ImSgn;


// Complex field
C<I> := ComplexField(20);
RPI := Real(Pi(C));


// Edges
declare type SEEdge;

declare attributes SEEdge :
	EP, // End points = [E1,E2]
	Data, // Data = [ u_1, ... , u_{n-2},(b-a)/2,(b+a)/(b-a),C_ab ]
	up,
	r,
	Vr,
	Isgn,
	IntSch,
	vp;

intrinsic Print( E::SEEdge )
{ Printing }
	print "Endpoints:",E`EP;
	/*print "Data:",E`Data;
	print "up:",E`up;
	print "r:",E`r;
	print "IntSch:",E`IntSch;*/
end intrinsic;

// Constructor
function SE_Edge(k,l:r:=0)
	E := New(SEEdge);
	if Type(l) eq RngIntElt then
		E`EP := [k,l];
	elif Type(l) eq SeqEnum then
		E`EP := <k,l>;
	else
		error Error("Not supposed to happen.");
	end if;
	if r ne 0 then
		E`r := r;
	end if;
	return E;
end function;

procedure EdgeData( E,Points,Zetas,m,n )
	a := Points[E`EP[1]]; 
	if Type(E`EP[2]) eq RngIntElt then
		b := Points[E`EP[2]];
		if E`EP[1] lt E`EP[2] then
			Pts := Remove(Remove(Points,E`EP[1]),E`EP[2]-1);
		else
			Pts := Remove(Remove(Points,E`EP[2]),E`EP[1]-1);
		end if;
		bpa := b+a; bma := b-a; bmainv := 1/bma;
		CCV, up := SortByRealPart([ (2*x-bpa)*bmainv : x in Pts ]);
		E`Isgn := [ ImSgn(CCV[k]) : k in [1..up] ] cat [ -ImSgn(CCV[k]) : k in [up+1..n-2] ];
		Append(~CCV,bma/2);
		Append(~CCV,(b+a)/bma);
		if IsReal(bma) and Real(bma) lt 0 then
			// Fix Magma here: Log(-R) inconsistent!
			C_ab := Exp( (n/m) * (Log(-bma)+Universe(Zetas).1*Pi(Parent(a)))) * Zetas[(up+1) mod 2 + 1];
		else
			C_ab := Exp( (n/m) * Log(bma)) * Zetas[(up+1) mod 2 + 1];
		end if;
		Append(~CCV,C_ab);
		E`up := up;
		E`Data := CCV;
	elif Type(E`EP[2]) eq SeqEnum then
		b := E`EP[2][1];
		Pts := Remove(Points,E`EP[1]);
		bpa := b+a; bma := b-a; bmainv := 1/bma;
		CCV, up := SortByRealPart([ (2*x-bpa)*bmainv : x in Pts ]);
		E`Isgn := [ ImSgn(CCV[k]) : k in [1..up] ] cat [ -ImSgn(CCV[k]) : k in [up+1..n-1] ];
		Append(~CCV,bma/2);
		Append(~CCV,(b+a)/bma);
		E`up := up;
		E`Data := CCV;
	else
		error Error("Not supposed to happen.");
	end if;
end procedure;

procedure FlipEdge( E )
	E`EP := Reverse(E`EP);
end procedure;

// Compare r-value of edges
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
procedure DE_AJM_Weight( ~Edge, Len : Lambda := RPI/2 )
	Edge`r := 5.0;
	for k in [1..Len] do
		Edge`r := Min(Edge`r,Abs(Im(Argsinh(Argtanh(Edge`Data[k])/Lambda))));
	end for;
	//Edge`r *:= 9/10;
end procedure;
procedure GJ_AJM_Weight( ~Edge, Len )
	Edge`Vr := [];
	for k in [1..Len] do
		P := ChangePrecision(Edge`Data[k],20);
		Append(~Edge`Vr,(Abs(P+1) + Abs(P-1))/2);
	end for;
	Edge`r := Min(Edge`Vr);
	//Append(~Edge,V_r);
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
	
	// ComplexEdges = <[P],v_P,k,r_E>
	m := SEC`Degree[1]; 
	n := SEC`Degree[2];
	NEdges := #ComplexEdges; 
	NewEdges := [];
	
	rm := 5.;
	for j in [1..NEdges] do
		NewEdge := SE_Edge(ComplexEdges[j][3],ComplexEdges[j][1]);
		NewEdge`vp := ComplexEdges[j][2];
		EdgeData(NewEdge,SEC`BranchPoints,[],m,n);
		GJ_AJM_Weight(~NewEdge,n-1);
		rm := Min(rm,NewEdge`r);
		Append(~NewEdges,NewEdge);
	end for;

	eps := 1/500;
	IntPars := [ NewEdges[k]`r : k in [1..NEdges] ];
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
		l := Max([ l : l in [1..NSchemes] | NewEdges[k]`r gt Lr[l] ]);
		NewEdges[k]`IntSch := l;
		M := Exp( -(1/m) * Log( &*[ NewEdges[k]`Vr[j] - Lr[l] : j in [1..SEC`Degree[2]-1]]));
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
	rmi := STree`IntPars[1];
	rma := STree`IntPars[2]; 
	vprint SE,2 : "Min_r:",rmi;
	vprint SE,2 : "Max_r:",rma;	
	NSchemes := Max(Ceiling(20*(rma-rmi)),1);
	IntPars := [ STree`Edges[k]`r : k in [1..STree`Length ] ];
	Groups := [ [] : j in [1..NSchemes] ];
	for r in IntPars do
		t := Max(Ceiling(20*(r-rmi)),1);
		Append(~Groups[t],r);
	end for;
	vprint SE,2 : "Groups:",Groups;
	Lr := [(29/30) * rmi ] cat [ (29/30) * Min(g) : g in Remove(Groups,1) | #g gt 0 ];
	NSchemes := #Lr;
	vprint SE,2 : "Number of schemes:",NSchemes;
	assert &and[rma lt RPI/2, rmi gt 0];

	// Find best integration scheme for each edge and compute bounds M1,M2
	NSchemes := #Lr; 
	LrM := [ <1,1,Lr[l]> : l in [1..NSchemes] ];
	for k in [1..STree`Length] do
		l := NSchemes;
		while STree`Edges[k]`r lt Lr[l] do
			l -:= 1;
		end while;
		STree`Edges[k]`IntSch := l;
		M1 := Ceiling(Bound_M1(STree`Edges[k]`Data,STree`Length-1,m));
		M2 := Ceiling(Bound_M2(STree`Edges[k]`Data,STree`Length-1,m,STree`Length-1,Lr[l]));
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
	
	rmi := 5.;
	rma := 0.;	
	for j in [1..NEdges] do
		NewEdge := SE_Edge(ComplexEdges[j][3],ComplexEdges[j][1]);
		NewEdge`vp := ComplexEdges[j][2];
		EdgeData(NewEdge,SEC`BranchPoints,[],m,n);
		DE_AJM_Weight(~NewEdge,n-1);
		rmi := Min(rmi,NewEdge`r);
		rma := Max(rma,NewEdge`r);
		Append(~NewEdges,NewEdge);
	end for;
	// NewEdges = <[P],v_P,k,r_E,CCV>

	// Make list of r's
	vprint SE,2 : "Min_r:",rmi;
	vprint SE,2 : "Max_r:",rma;	
	NSchemes := Max(Ceiling(20*(rma-rmi)),1);
	IntPars := [ NewEdges[k]`r : k in [1..NEdges ] ];
	Groups := [ [] : j in [1..NSchemes] ];
	for r in IntPars do
		t := Max(Ceiling(20*(r-rmi)),1);
		Append(~Groups[t],r);
	end for;
	vprint SE,2 : "Groups:",Groups;
	Lr := [ (29/30) * rmi ] cat [ (29/30) * Min(g) : g in Remove(Groups,1) | #g gt 0 ];
	NSchemes := #Lr;
	vprint SE,2 : "Number of schemes:",NSchemes;
	assert &and[rma lt RPI/2, rmi gt 0];


	// Find best integration scheme for each edge
	NSchemes := #Lr; 
	LrM := [ <1,1,Lr[l]> : l in [1..NSchemes] ];
	for k in [1..NEdges] do
		l := NSchemes;
		while NewEdges[k]`r lt Lr[l] do
			l -:= 1;
		end while;
		NewEdges[k]`IntSch := l;
		M1 := Ceiling(Bound_M1(NewEdges[k]`Data,n-1,m:AJM));
		M2 := Ceiling(Bound_M2(NewEdges[k]`Data,n-1,m,n-2,Lr[l]:AJM));
		LrM[l][1] := Max(M1,LrM[l][1]);
		LrM[l][2] := Max(M2,LrM[l][2]);
	end for;

	vprint SE,2 : "Parameters(AJM):",LrM;

	// Return parameters
	return LrM, NewEdges;
end function;


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



