/*******************************************************************************

 Copyright (C) 2017 
 Code adapted from Pascal Molin,
 written by Christian Neurohr

*******************************************************************************/

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
	if Type(E[2]) eq RngIntElt then
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
// Computes the roots of the complex polynomial f to D decimal digits from approximations Z
	f *:= (1/LeadingCoefficient(f));
	f := ChangeRing(f,ComplexField(2*Digits));
	N := Degree(f);
	RMV := [ Remove([1..N],j) : j in [1..N] ];
	Err2 := (1/2) * 10^-(Digits+1);
	W := [ Evaluate(f,Z[j])/ &*[ (Z[j] - Z[k]) : k in RMV[j] ] : j in [1..N] ];
	w0 :=  Max([ Abs(W[j]) : j in [1..N] ]);
	if w0 lt Err2 then
		return ChangeUniverse(Z,ComplexField(Digits));
	end if;
	p := Precision(Universe(Z));
	d0 := DistanceII(Z);
	if 2*w0 lt d0 then
		repeat
			Z := [ Z[j] - W[j] : j in [1..N] ];
			p := Max(2*p,Digits);
			ChangeUniverse(~Z,ComplexField(p));
			W := [ Evaluate(f,Z[j])/ &*[ (Z[j] - Z[k]) : k in RMV[j] ] : j in [1..N] ];
			w0 := Max([ Abs(W[j]) : j in [1..N] ]);
		until w0 lt Err2;
		return Sort(Z,SE_OrdFldComElt);
	else
		assert false;
		return [];
	end if;
end function;





