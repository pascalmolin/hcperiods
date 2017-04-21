/******************************************************************************

 Copyright (C) 2017 
 Code adapted from William Stein,
 written by Christian Neurohr

 ******************************************************************************/

procedure RS_InplaceMoveToPositivePivot( Column, Row, Pivot, ~G, ~B )
	i := Column; j := Row; v := G[i][j]; IsPivot := false;
	if [i,j] eq [Pivot+1,Pivot] and G[Pivot+1][Pivot] ne v then
		IsPivot := true;
        	SwapRows(~B,Pivot,Pivot+1);
        	SwapRows(~G,Pivot,Pivot+1);
        	SwapColumns(~G,Pivot,Pivot+1);
    	elif [i,j] eq [Pivot,Pivot+1] then
        	SwapRows(~B,Pivot,Pivot+1);
        	SwapRows(~G,Pivot,Pivot+1);
        	SwapColumns(~G,Pivot,Pivot+1);
    	elif j ne Pivot and j ne Pivot+1 and i ne Pivot and i ne Pivot+1 then
        	SwapRows(~B,Pivot,j);
	        SwapRows(~B,Pivot+1,i);
		SwapRows(~G,Pivot,j);
        	SwapRows(~G,Pivot+1,i);
        	SwapColumns(~G,Pivot,j);
        	SwapColumns(~G,Pivot+1,i);
    	elif j eq Pivot then
        	SwapRows(~B,Pivot+1,i);
        	SwapRows(~G,Pivot+1,i);
        	SwapColumns(~G,Pivot+1,i);
   	elif j eq Pivot+1 then
        	SwapRows(~B,Pivot,i);
        	SwapRows(~G,Pivot,i);
        	SwapColumns(~G,Pivot,i);
    	elif i eq Pivot then
        	SwapRows(~B,Pivot+1,j);
        	SwapRows(~G,Pivot+1,j);
        	SwapColumns(~G,Pivot+1,j);
    	elif i eq Pivot+1 then
        	SwapRows(~B,Pivot,j);
      		SwapRows(~G,Pivot,j);
        	SwapColumns(~G,Pivot,j);
	end if;
    	// Fix possible change of sign in a row
    	if G[Pivot+1][Pivot] ne v and not IsPivot then
        	SwapRows(~B,Pivot,Pivot+1);
        	SwapRows(~G,Pivot,Pivot+1);
        	SwapColumns(~G,Pivot,Pivot+1);
	end if;
end procedure;


function RS_FindSmallestElementPosition( K, Pivot )
	n := NumberOfRows(K);
	for i in [Pivot..n] do
		for j in [Pivot..n] do
			if K[i][j] eq 1 then
				return [i,j];
			end if;
		end for;
	end for;
	return [0,0];
end function;


function SymplecticBasis(K)
	E := K;
	B := One(Parent(E));
    	Zeros := [];
    	PS := [];
    	Pivot := 1;
	n := NumberOfColumns(K);
    	while Pivot le n do
        	SmallestEltPos := RS_FindSmallestElementPosition( E, Pivot );
        	if SmallestEltPos eq [0,0] then
            		Append(~Zeros,Pivot);
           		Pivot := Pivot + 1;
            		continue;
		end if;
       		RS_InplaceMoveToPositivePivot(SmallestEltPos[2], SmallestEltPos[1], Pivot, ~E, ~B);
		AllZero := true;

		// Use non-zero element to clean row Pivot
       	       	u := E[Pivot][Pivot+1];
       		for j in [Pivot+2..n] do
        		v, r := Quotrem(-E[Pivot][j],u);
            		if v ne 0 then
                		AllZero := false;
                		AddRow(~E,v,Pivot+1,j);
                		AddColumn(~E,v,Pivot+1,j);
               			AddRow(~B,v,Pivot+1,j);
			end if;
		end for;
		
        	// Use non-zero element to clean row Pivot+1
        	u := E[Pivot+1][Pivot];	
	     	for j in [Pivot+2..n] do
        		v, r := Quotrem(-E[Pivot+1][j],u);
            		if v ne 0 then
               			AllZero := false;
				AddRow(~E,v,Pivot,j);
				AddColumn(~E,v,Pivot,j);
				AddRow(~B,v,Pivot,j);
			end if;
		end for;
		
		// Record for basis reconstruction
           	if AllZero then
            		Append(~PS,[E[Pivot+1][Pivot], Pivot]);
            		Pivot := Pivot +  2;
		end if;
	end while;
	
    	Sort(~PS);
	Reverse(~PS);
	
  	ES := [ p[2] : p in PS ];
    	FS := [ p[2]+1 : p in PS ];
	
	NewRowsIndices := ES cat FS cat Zeros;
	
	NewRows := [ RowSubmatrix(B,i,1) : i in NewRowsIndices ];
    	ST := Matrix(NewRows); 		// Base change matrix
  	CF := ST * K * Transpose(ST);   // Symplectic basis

    	return CF, ST;
end function;
