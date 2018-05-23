/*******************************************************************************

 Copyright (C) 2017
 Christian Neurohr

*******************************************************************************/

function ImSgn(z)
	if IsReal(z) then
		return 1;
	else
		return Sign(Im(z));
	end if;
end function;

function AC_mthRoot(x,E,Zetas,m,n)
// Used for analytic continuation of y = f(x)^(1/m)
	WindingNr := 0;
	if 1 le E`up then
		z := E`Data[1]-x;
	else
		z := x-E`Data[1];
	end if;
	isgnz := E`Isgn[1];
	for k in [2..n] do
		if k le E`up then
			z *:= E`Data[k]-x;
		else
			z *:= x-E`Data[k];
		end if;
		if isgnz eq E`Isgn[k] then
			isgnz := ImSgn(z);
			if isgnz ne E`Isgn[k] then
				if E`Isgn[k] eq 1 then
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
		return Zetas[WindingNr mod 4 + 1]  * Sqrt(z);
	else
		return Zetas[WindingNr mod (2*m) + 1]  * Exp(Log(z)/m);
	end if;
end function;
function AC_mthRoot2(x,E,Zetas,m,n)
// Used for analytic continuation of y = f(x)^(1/m)
	z := 0;
	for k in [1..E`up] do
		z +:= Log(E`Data[k] - x);
	end for;
	for k in [E`up+1..n] do
		z +:= Log(x - E`Data[k]);
	end for;
	return Exp(z/m);
end function;
function ACRec(dat,sgns,D)
		//D := #dat;
		//assert IsCoercible(Integers(),Log(2,D));
		/*if D eq 1 then
			return dat[1],sgns[1], w;*/
		//assert D gt 1;
		w := 0;
		if D eq 2 then
			z := dat[1]*dat[2];
			if sgns[1] eq sgns[2] then
				isgnz := ImSgn(z);
				if isgnz ne sgns[2] then
					if sgns[2] eq 1 then
						w +:= 2;
					else
						w -:= 2;
					end if;
				end if;
			else
				isgnz := ImSgn(z);
			end if;
			return z, isgnz, w;
		else
			D2 := Round(D/2);//print "sgns:",sgns;print "D2:",D2;print "D:",D;print "w:",w;
			z1,sgn1,w1 := ACRec(dat[1..D2],sgns[1..D2],D2);
			z2,sgn2,w2 := ACRec(dat[D2+1..D],sgns[D2+1..D],D2);
			w := w1+w2;
			z := z1*z2;
			if sgn1 eq sgn2 then
				isgnz := ImSgn(z);
				if isgnz ne sgn2 then
					if sgn2 eq 1 then
						w +:= 2;
					else
						w -:= 2;
					end if;
				end if;
			else
				isgnz := ImSgn(z);
			end if;
			return z, isgnz, w;
		end if;
	end function;
function AC_mthRoot3(x,E,Zetas,m,n)
// Used for analytic continuation of y = f(x)^(1/m)
	w := 0;
 	l := Ceiling(Log(2,n));
	Data := [ d-x : d in E`Data[1..E`up] ] cat [ x-d : d in E`Data[E`up+1..n] ] cat [ 1 : j in [n+1..2^l] ];
	Sgns := [ e : e in E`Isgn[1..E`up] ] cat [ -e : e in E`Isgn[E`up+1..n] ] cat [ 1 : j in [n+1..2^l] ];
	D := #Data;
	if D gt 1 then
		z, sgnz, w := ACRec(Data,Sgns,D);
	else
		z := Data[1];
	end if;
	if Im(z) eq 0 and Re(z) lt 0 then
		z := -z;
		w +:= 1;	
	end if;
	if m eq 2 then
		return Zetas[w mod 4 + 1]  * Sqrt(z);
	else
		return Zetas[w mod (2*m) + 1]  * Exp(Log(z)/m);
	end if;
end function;

