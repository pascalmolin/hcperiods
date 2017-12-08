/*******************************************************************************

 Copyright (C) 2017 
 Code adapted from Pascal Molin,
 written by Christian Neurohr

*******************************************************************************/

function IntGroups(IntParPaths,rm,m)
	if m gt 2 then
		Groups := [ [],[],[],[] ];
		for r in IntParPaths do
			if r lt rm+0.4 then Append(~Groups[1],r);
			else Append(~Groups[4],r); end if;
		end for;
	else
		Groups := [ [],[],[],[],[],[] ];
		for r in IntParPaths do if r lt rm+0.1 then Append(~Groups[1],r); elif r lt rm+0.2 then Append(~Groups[2],r);
			elif r lt rm+0.4 then Append(~Groups[3],r);
			elif r lt rm+0.8 then Append(~Groups[4],r);
			elif r lt rm+1.6 then Append(~Groups[5],r);
			else Append(~Groups[6],r); end if;
		end for;
	end if;
	return Groups;
end function;

function SortByRealPart(V)
	oV := []; up := 0;
	for z in V do
		if Re(z) le 0 then
			Append(~oV,z);	
		else
			Insert(~oV,1,z);
			up +:= 1;
		end if;
	end for;
	return oV,up;
end function;

function MakeCCVector( E, Points )
	// E = < a, b > or E = < i , b >
	assert Type(E[1]) eq RngIntElt;
	a := Points[E[1]];
	if Type(E) eq SEEdge then
		b := Points[E[2]];
		if E[1] lt E[2] then
			Pts := Remove(Remove(Points,E[1]),E[2]-1);
		else
			Pts := Remove(Remove(Points,E[2]),E[1]-1);
		end if;
		CCV, up := SortByRealPart([ (2*x-b-a)/(b-a) : x in Pts ]);
		bma := b-a;
		Append(~CCV,bma/2);
		Append(~CCV,(b+a)/bma);
	elif Type(E[2]) eq FldComElt then
		a := Points[E[1]];
		p := E[2];
		Pts := Remove(Points,E[1]);
		CCV, up := SortByRealPart([ (2*x-p-a)/(p-a) : x in Pts ]);
		pma := p-a;
		Append(~CCV,pma/2);
		Append(~CCV,(p+a)/pma);
	else
		error Error("Not supposed to happen.");
	end if;
	return CCV, up;
end function;


procedure PolynomialShiftVector( ~V, c, Len, Shft )
// returns (v0, c*v0 + v1, c^2*v0 + 2c*v1 + v2, ...)
	for k in [2..Len] do
		l := Len;
		while l ge k do
			V[l+Shft] +:= c * V[l+Shft-1];
			l -:= 1;
		end while;
	end for;
end procedure;


function Distance(x,Points)
	return Min([ Abs(x-P) : P in Points ]);
end function;
function DistanceII(Points)
	Distances := [];
	L := #Points;
	for k in [1..L] do
		for l in [k+1..L] do
			Append(~Distances,Abs(Points[k]-Points[l]));
		end for;
	end for;
	return Min(Distances);
end function;

function SE_OrdFldComElt(x,y)
// Used to fix an ordering of the sheets
	if Re(x) lt Re(y) then 
		return -1;
	elif Re(x) gt Re(y) then 
		return 1; 
	else 
		if Im(x) lt Im(y) then 
			return -1;
		elif Im(x) gt Im(y) then 
			return 1; 
		else
			return 0;
		end if; 
	end if;
end function;


function SE_DKPEB( f,Z,Digits )
	f := ChangeRing(f/LeadingCoefficient(f),ComplexField(2*Digits));
	Z := ChangeUniverse(Z,ComplexField(2*Digits));
	N := Degree(f);
	RMV := [ Remove([1..N],j) : j in [1..N] ];
	Err2 := (1/2) * 10^-(Digits+1);
	// Start root approximation
	W := [ Evaluate(f,Z[j])/ &*[ (Z[j] - Z[k]) : k in RMV[j] ] : j in [1..N] ];
	repeat
		Z := [ Z[j] - W[j] : j in [1..N] ];
		W := [ Evaluate(f,Z[j])/ &*[ (Z[j] - Z[k]) : k in RMV[j] ] : j in [1..N] ];
		w0 := Max([ Abs(W[j]) : j in [1..N] ]);
	until w0 lt Err2;
	return Z;
end function;




