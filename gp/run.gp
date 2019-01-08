install("gc_init","Lp","gc_init","./hyperellperiods.so");
install("hcinit","Gp","hcinit","./hyperellperiods.so");
install("symplectic_reduction","GLD&","symplectic_reduction","./hyperellperiods.so");
{
/* chebychev nodes */
for(n=2,30,
  gc = gc_init(n);
  ref = [ cos((2*k-1)*Pi/(4*n)) | k <- [1..n] ];
  if(exponent(vecmax(abs(gc-ref)))>-110,error("rootsof1 ",n,"\n",gc,"\n---\n",ref));
  );
/* hcinit */
pols = [x^4-1,polcyclo(13),x^5+2*x+3,polhermite(11),polcyclo(18)];
for(i=1,#pols,
  hc = hcinit(pols[i]);
  );
}
{
/* symplectic reduction */
for(g=1,20,
  for(n=1,5,
    my(m=matrix(2*g,2*g));
    for(i=1,2*g,for(j=i+1,2*g,if(random(2),m[i,j]=-m[j,i]=random(5)-5)));
    my(p=symplectic_reduction(m,g,&d));
    jg=matconcat(matdiagonal(vector(g,k,d[k]*[0,1;-1,0])));
    if(p~*m*p != jg, error(m));
    ));
}
