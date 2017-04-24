/*******************************************************************************

 Copyright (C) 2017
 Christian Neurohr

*******************************************************************************/


function ImSgn( z )
	if IsReal(z) then
		return 1;
	else
		return Sign(Im(z));
	end if;
end function;
function AC_mthRoot(x,CCV,Zetas,up,m,n)
// Used for analytic continuation of y = f(x)^(1/m)
	WindingNr := 0;
	if 1 le up then
		z := CCV[1] - x;
	else
		z := x - CCV[1];
	end if;
	isgnz := ImSgn(z);
	for k in [2..n] do
		if k le up then
			z_ := CCV[k] - x;
		else
			z_ := x - CCV[k];
		end if;
		isgnz_ := ImSgn(z_);
		z *:= z_;
		if isgnz eq isgnz_ then
			isgnz := ImSgn(z);
			if isgnz ne isgnz_ then
				if isgnz_ gt 0 then
					WindingNr +:= 2;
				else
					WindingNr -:= 2;
				end if;
			end if;
		else
			isgnz := ImSgn(z);
		end if;
	end for;
	if Im(z) eq 0 and Re(z) lt 0 then
		z := -z;
		WindingNr +:= 1;	
	end if;
	if m eq 2 then
		z := Zetas[WindingNr mod (2*m) + 1]  * Sqrt(z);
	else
		z := Zetas[WindingNr mod (2*m) + 1]  * Exp((1/m)*Log(z));
	end if;
	return z;
end function;
