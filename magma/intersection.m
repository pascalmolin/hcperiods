//Import global settings
import "globalsettings.m": RS_Config;
import "superelliptic.m": RS_ImSgn, RS_NRootAC;
import "comparefunctions.m": RS_CompareFldComElt;


intrinsic RS_SEIntersection( Edge1::Tup, Edge2::Tup, Points::[FldComElt], N::RngIntElt ) -> Tup
{ - }
	NrOfEndPts := #Set([Edge1[1],Edge1[2],Edge2[1],Edge2[2]]);
	if NrOfEndPts eq 4 then
		return <0,0,0>; // cannot occur in MST
		//print "################################ Case5: nce ################################";
		//return RS_InnerIntersection( Edge1,Edge2,Points,N );
	end if;

	// Trivial cases
	if Edge1[1] eq Edge1[2] or Edge2[1] eq Edge2[2] then
		// a = b or c = d
		error Error("Bad edge in MST.");
	end if;
	if ( Edge1[1] eq Edge2[1] and Edge1[2] eq Edge2[2] ) or ( Edge1[1] eq Edge2[2] and Edge1[2] eq Edge2[1] ) then
		// ( a = c and b = d ) or ( a = d and b = c )
		//not true for N > 2// will be replaced later
		"Case0: Self-intersection";
		error Error("How can this happen?.");
		return <0,0,0>;
	end if;

	// Number of points
	l := #Points;
	assert l ge 3;

	// Variables;
	a := Points[Edge1[1]]; b := Points[Edge1[2]];
	c := Points[Edge2[1]]; d := Points[Edge2[2]];	
	Zeta := RS_GetGlobalZeta(); Pi := RS_GetGlobalPi();

	// Constants
	c1 := (b-a)^(l/N); 
	c2 := (d-c)^(l/N);

	// Angle between linesegments
	phi := Arg((d-c)/(b-a));
	

	assert &and[-Pi lt phi+10^(-10), phi-10^(-10) lt Pi];
	//assert &and[-Pi lt phi, phi lt Pi];

	// Making vectors of centers of circumcircles
	PPP1 := []; PPP2 := [];
	for k in [1..l] do
		if k notin [Edge1[1],Edge1[2]] then
			Append(~PPP1,(2*Points[k]-b-a)/(b-a));
		end if;
		if k notin [Edge2[1],Edge2[2]] then
			Append(~PPP2,(2*Points[k]-d-c)/(d-c));
		end if;
	end for;

	/*
	print "########################### RS_SEIntsec ############################";
	print "Edge1:",Edge1; print "Edge2:",Edge2;
	print "phi:",phi; 
	print "Points:",Points;
	print "#Points:",#Points;
	print "N:",N;
	print "c1:",c1;
	print "c2:",c2;
	print "phi:",phi; 
	print "PPP1:",PPP1;
	print "PPP2:",PPP2;
	*/

	// Cases
	if Edge1[1] eq Edge2[1] then
		//"################################ Case1: a = c ################################";
		AR1 := RS_NRootAC(-1,PPP1,Zeta,N);
		AR2 := RS_NRootAC(-1,PPP2,Zeta,N);
		Val_ab := c1 * AR1;
		Val_ad := c2 * AR2;		
		Val := Val_ab/Val_ad;
		ValArg := Arg(Val);
		//return <Sign(Im(Val)),0>;
		
		if Abs(phi) le Pi/2 then
			LimArg := -phi/2;
			ints := -Sign(phi);
			if phi gt 0 then
				dir := 1;
			else
				dir := N-1;
			end if;
		elif phi gt 0 then
			LimArg := Pi - phi/2;
			ints := 1;
			dir := N-1;
		elif phi lt 0 then
			LimArg := -Pi - phi/2;
			ints := -1;
			dir := 1;
		end if;
		//ints := -1;
		
	elif Edge1[2] eq Edge2[1] then
		//"################################ Case2: b = c ################################";
		AR1 := RS_NRootAC(1,PPP1,Zeta,N);
		AR2 := RS_NRootAC(-1,PPP2,Zeta,N);
		Val_ab := c1 * AR1;
		Val_ad := c2 * AR2;
		/*
		Val := Val_ad/Val_ab;
		ValArg := Arg(Val);
		print "ValArg:",ValArg;
		return <Sign(Im(Val)),0>;
		*/
		Val := Val_ab/Val_ad;
		ValArg := Arg(Val);
		
		if phi lt 0 then
			LimArg := (Pi-phi)/2;
			ints := -1;
		else
			LimArg := -(Pi+phi)/2;
			ints := 1;
		end if;
		
		LimArg := (Pi-phi)/2;
		ints := -1;
		dir := N-1;
	elif Edge1[1] eq Edge2[2] then
		//"################################ Case3: a = d ################################";
		//sk := RS_SEIntersection( Edge2, Edge1, Points, N );
		//return <-sk[1],sk[2]>;
		AR1 := RS_NRootAC(-1,PPP1,Zeta,N);
		AR2 := RS_NRootAC(1,PPP2,Zeta,N);
		Val_ab := c1 * AR1;
		Val_cd := c2 * AR2;		
		Val := Val_ab/Val_cd;
		ValArg := Arg(Val);
		LimArg := (Pi-phi)/2;
		ints := -1;
		dir := N-1;
	elif Edge1[2] eq Edge2[2] then
		//"################################ Case1: b = d ################################";
		AR1 := RS_NRootAC(1,PPP1,Zeta,N);
		AR2 := RS_NRootAC(1,PPP2,Zeta,N);
		Val_ab := c1 * AR1;
		Val_cd := c2 * AR2;		
		Val := Val_ab/Val_cd;
		//Val := Val_cd/Val_ab;
		ValArg := Arg(Val);
		//return <Sign(Im(Val)),0>;
		
		if Abs(phi) le Pi/2 then
			LimArg := -phi/2;
			ints := -Sign(phi);
			if phi gt 0 then
				dir := 1;
			else
				dir := N-1;
			end if;
		elif phi gt 0 then
			LimArg := Pi - phi/2;
			ints := 1;
			dir := N-1;
		elif phi lt 0 then
			LimArg := -Pi - phi/2;
			ints := -1;
			dir := 1;
		end if;
	else
		//"################################ Case5: nce ################################";
		return <0,0,0>; // cannot occur in MST
		//return RS_InnerIntersection( Edge1,Edge2,Points,N );
	end if;
	
	// Compute shift between sheets
	k := (1/Pi) * ( LimArg - (N/2) * ValArg );
	k_ := Round(k);
	assert Abs(k - k_) lt RS_GetGlobalError(); // should be close to an integer
	k := k_ mod N;
	return <ints,k,dir>;

	
end intrinsic;



intrinsic RS_InnerIntersection( Edge1::Tup,Edge2::Tup,Points::[FldComElt],N::RngIntElt ) -> RngIntElt
{ - }
	C<i> := Parent(Points[1]); Pi := RS_GetGlobalPi(); RS_SetGlobalZeta(N); Zeta := RS_GetGlobalZeta(); l := #Points;
	a := Points[Edge1[1]]; b := Points[Edge1[2]];
	c := Points[Edge2[1]]; d := Points[Edge2[2]];

	// Check if Int(Edge1,Edge2) = 0
	cprim := (2*c - a - b) / (b-a);
	dprim := (2*d - a - b) / (b-a);
	imc := Im(cprim);
	imd := Im(dprim);
	if Sign(imc)*Sign(imd) eq 1 then
		//print "// on the same side";
		"!!!!!!!!!!!!!!!!!!!!";
		return 0; 
	end if;
	xp := Im(Conjugate(cprim)*(dprim-cprim));
	fp := imd - imc;
	/*
	print "#Points:",#Points;
	print "Edge1 from",a,"to",b;
	print "Edge2 from",c,"to",d;
	print "xp:",xp;
	print "fp:",fp;
	*/
	print "Sign(fp) (IntNr of [ab].[cd]):",Sign(fp);

	if Abs(xp) ge Abs(fp) then // discard, if xp not in ]-1,1[
		"!!!!!!!!!!!!!!!!!!!!!";
		return 0;
	end if;

	xpab := xp/fp;
	p := (b+a)/2 + xpab*((b-a)/2);
	xpcd := (2*p - c - d)/(d-c); assert &and[-1 lt Real(xpcd),Real(xpcd) lt 1];
	print "Intersection point:",p;

	PPP1 := []; PPP2 := [];
	for k in [1..#Points] do
		if k notin [Edge1[1],Edge1[2]] then
			Append(~PPP1,(2*Points[k]-b-a)/(b-a));
		end if;
		if k notin [Edge2[1],Edge2[2]] then
			Append(~PPP2,(2*Points[k]-d-c)/(d-c));
		end if;
	end for;
	x_alpha_p := (2*p - a - b)/(b-a); assert &and[-1 lt Real(x_alpha_p),Real(x_alpha_p) lt 1];
	w_alpha_p := Argtanh(x_alpha_p); 
	x_beta_p := (2*p - c - d)/(d-c); assert &and[-1 lt Real(x_beta_p),Real(x_beta_p) lt 1];
	w_beta_p := Argtanh(x_beta_p); 
	AR1 := RS_NRootAC(x_alpha_p,PPP1,Zeta,N); 
	AR2 := RS_NRootAC(x_beta_p,PPP2,Zeta,N);
	assert &and[Real(Cosh(w_alpha_p)) gt 0, Real(Cosh(w_beta_p)) gt 0];
	f_alpha_p := (Cosh(w_alpha_p)^2 * Exp(-Pi*i) * (2/(b-a))^l)^(1/N) / AR1;
	f_beta_p := (Cosh(w_beta_p)^2 * Exp(-Pi*i) * (2/(d-c))^l)^(1/N) / AR2;
	f_alpha_p := (Cosh(w_alpha_p))^(2/N) * Exp(Pi*i/N) * (2/(b-a))^(l/N) / AR1;
	f_beta_p := (Cosh(w_beta_p))^(2/N) * Exp(Pi*i/N) * (2/(d-c))^(l/N) / AR2;
	f_alpha_p := (Cosh(w_alpha_p))^(2/N) * (2/(b-a))^(l/N) / AR1;
	f_beta_p := (Cosh(w_beta_p))^(2/N)  * (2/(d-c))^(l/N) / AR2;
	quot := f_alpha_p/f_beta_p;
	quot := (Cosh(w_alpha_p)/Cosh(w_beta_p))^(2/N) * ((d-c)/(b-a))^(l/N) * (AR2/AR1);
	print "quot(RoU):",quot;
	print "Abs(quot):",Abs(quot);
	/*
	xi := (Exp(-Pi*i) * (2/(b-a))^l)^(1/N) / (Exp(-Pi*i) * (2/(d-c))^l)^(1/N) /  ((d-c)/(b-a))^(l/N);
	print "xi:",xi;
	quot2 := xi * (Cosh(w_alpha_p)/Cosh(w_beta_p))^(2/N) * ((d-c)/(b-a))^(l/N) * (AR2/AR1);
	print "quot2:",quot2;
	*/
	k_ := (N/(2*Pi)) * Arg(quot);
	print "k_:",k_;
	k := Round(k_) mod N;
	print "k:",k;
	return k;

	return 2*Sign(fp)*Sign(Real(f_alpha_p/f_beta_p));
end intrinsic;



intrinsic RS_IntervalIntersection( Edge1::Tup,Edge2::Tup,Points::[FldComElt] ) -> RngIntElt
{ - }
	C<i> := Parent(Points[1]); Pi := RS_GetGlobalPi(); l := #Points;
	a := Points[Edge1[1]]; b := Points[Edge1[2]];
	c := Points[Edge2[1]]; d := Points[Edge2[2]];

	// Check if Int(Edge1,Edge2) = 0
	cprim := (2*c - a - b) / (b-a);
	dprim := (2*d - a - b) / (b-a);
	imc := Im(cprim);
	imd := Im(dprim);
	if Sign(imc)*Sign(imd) eq 1 then
		//print "// on the same side";
		return 0; 
	end if;
	xp := Im(Conjugate(cprim)*(dprim-cprim));
	fp := imd - imc;
	/*
	print "#Points:",#Points;
	print "a:",a;
	print "b:",b;
	print "c:",c;
	print "d:",d;
	print "fp:",fp;
	print "Sign(fp) (IntNr of [ab].[cd]):",Sign(fp);
	*/
	if Abs(xp) ge Abs(fp) then // discard, if xp not in ]-1,1[
		return 0;
	end if;
	
	return Sign(fp);
end intrinsic;


