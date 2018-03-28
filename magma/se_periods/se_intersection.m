/*******************************************************************************

 Copyright (C) 2017 
 Christian Neurohr

*******************************************************************************/

// Import functions
import "se_anal_cont.m": AC_mthRoot;
import "se_help_funcs.m": MakeCCVector;

PI302INV := 1/(2*Pi(RealField(30)));

function SE_IntersectionMatrix( spsm_Matrix, m,n )
	// Building block matrices
	Blocks := [];
	// Block matrix for self-shifts build the diagonal of intersection matrix
	SelfShiftBlock := ZeroMatrix(Integers(),m-1,m-1);
	for Ind1 in [1..m-2] do
		SelfShiftBlock[Ind1][Ind1+1] := 1;
		SelfShiftBlock[Ind1+1][Ind1] := -1;
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
	NrOfEndPts := #Set(Edge1`EP cat Edge2`EP);
	if NrOfEndPts eq 4 then
		return <0,0>;
	end if;

	// End points
	a := Points[Edge1`EP[1]]; b := Points[Edge1`EP[2]];
	c := Points[Edge2`EP[1]]; d := Points[Edge2`EP[2]];	

	// Trivial cases
	if a eq b or c eq d then
		error Error("Bad edge in spanning tree.");
	end if;
	if ( a eq c and b eq d ) or ( a eq d and d eq c ) then
		error Error("Bad edge in spanning tree.");
	end if;

	// Arg of factor corresponding to end points
	phi := Arg(Edge1`Data[n-1]/Edge2`Data[n-1]);

	// Cases
	if a eq c then
		//"################################ Case1: a = c ################################";
		AR1 := AC_mthRoot(-1,Edge1,Zetas,m,n-2);
		AR2 := AC_mthRoot(-1,Edge2,Zetas,m,n-2);
		Val_ab := Edge1`Data[n+1] * AR1;
		Val_cd := Edge2`Data[n+1] * AR2;		
		Val := Val_cd/Val_ab;
		k := Round( PI302INV * ( phi + m * Arg(Val) ));
		/*k := Round(k_);
		assert Abs(k-k_) lt 10^-10;*/
		if phi gt 0 then
			return <1-k,-k>;
		else
			return <-k,-1-k>;
		end if;
	elif b eq c then
		//"################################ Case2: b = c ################################";
		AR1 := AC_mthRoot(1,Edge1,Zetas,m,n-2);
		AR2 := AC_mthRoot(-1,Edge2,Zetas,m,n-2);
		Val_ab := Edge1`Data[n+1] * AR1;
		Val_cd := Edge2`Data[n+1] * AR2;		
		Val := Val_cd/Val_ab;
		k := Round( PI302INV * ( phi + m * Arg(Val) ) + (1/2));
		/*k := Round(k_);
		assert Abs(k-k_) lt 10^-10;*/
		return <-k,1-k>;
	else
		error Error("Bad edge in spanning tree.");
	end if;
end function;
