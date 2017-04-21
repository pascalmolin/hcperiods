/******************************************************************************

 Copyright (C) 2017 Christian Neurohr

 ******************************************************************************/

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
	// E = < a, b >
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









