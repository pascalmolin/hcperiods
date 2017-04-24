/*******************************************************************************

 Copyright (C) 2017 
 Christian Neurohr

*******************************************************************************/

// Import functions
import "se_anal_cont.m": AC_mthRoot;
import "se_help_funcs.m": MakeCCVector;


Pi30 := Pi(RealField(30));

function SE_IntersectionMatrix( spsm_Matrix, m,n )
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
	return BlockMatrix(n-1,n-1,Blocks);
end function;


function SE_IntersectionNumber( Edge1,Edge2,Points,m,n,Zetas ) 
	NrOfEndPts := #Set([Edge1[1],Edge1[2],Edge2[1],Edge2[2]]);
	if NrOfEndPts eq 4 then
		return <0,0>;
	end if;

	// Trivial cases
	if Edge1[1] eq Edge1[2] or Edge2[1] eq Edge2[2] then
		// a = b or c = d
		error Error("Bad edge in spanning tree.");
	end if;
	if ( Edge1[1] eq Edge2[1] and Edge1[2] eq Edge2[2] ) or ( Edge1[1] eq Edge2[2] and Edge1[2] eq Edge2[1] ) then
		// ( a = c and b = d ) or ( a = d and b = c )
		error Error("Bad edge in spanning tree.");
		return <0,0,0>;
	end if;

	// End points
	a := Points[Edge1[1]]; b := Points[Edge1[2]];
	c := Points[Edge2[1]]; d := Points[Edge2[2]];	

	// Arg of factor corresponding to end points
	phi := Arg((b-a)/(d-c));	

	// Making vectors of centers of circumcircles
	CCV1, up1 := MakeCCVector(Edge1,Points);
	CCV2, up2 := MakeCCVector(Edge2,Points);

	// Constants
	C_ab := Exp( (n/m) * Log(b-a)) * Zetas[(up1+1) mod 2 + 1]; 
	C_cd := Exp( (n/m) * Log(d-c)) * Zetas[(up2+1) mod 2 + 1];

	// Cases
	if Edge1[1] eq Edge2[1] then
		//"################################ Case1: a = c ################################";
		AR1 := AC_mthRoot(-1,CCV1,Zetas,up1,m,n-2);
		AR2 := AC_mthRoot(-1,CCV2,Zetas,up2,m,n-2);
		Val_ab := C_ab * AR1;
		Val_cd := C_cd * AR2;		
		Val := Val_cd/Val_ab;
		k_ := ((1/(2*Pi30)) * ( phi + m * Arg(Val) ));
		k := Round(k_);
		assert Abs(k-k_) lt 10^-10;
		if phi gt 0 then
			return <1-k,-k>;
		else
			return <-k,-1-k>;
		end if;
	elif Edge1[2] eq Edge2[1] then
		//"################################ Case2: b = c ################################";
		AR1 := AC_mthRoot(1,CCV1,Zetas,up1,m,n-2);
		AR2 := AC_mthRoot(-1,CCV2,Zetas,up2,m,n-2);
		Val_ab := C_ab * AR1;
		Val_cd := C_cd * AR2;		
		Val := Val_cd/Val_ab;
		k_ := (1/(2*Pi30)) * ( phi + m * Arg(Val) ) + (1/2);
		k := Round(k_);
		assert Abs(k-k_) lt 10^-10;
		return <-k,1-k>;
	else
		error Error("Bad edge in spanning tree.");
	end if;
end function;
