install("gc_init","Lp","gc_init","./hyperellperiods.so");
install("hcinit","Gp","hcinit","./hyperellperiods.so");
{
for(n=2,30,
  gc = gc_init(n);
  ref = [ cos((2*k-1)*Pi/(4*n)) | k <- [1..n] ];
  if(exponent(vecmax(abs(gc-ref)))>-110,error("rootsof1 ",n,"\n",gc,"\n---\n",ref));
  );
pols = [x^4-1,polcyclo(13),x^5+2*x+3,polhermite(11),polcyclo(18)];
for(i=1,#pols,
  hc = hcinit(pols[i]);
  );
}
