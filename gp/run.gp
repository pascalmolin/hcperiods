install("gc_init","Lp","gc_init","./hyperellperiods.so");
install("hcinit","Gp","hcinit","./hyperellperiods.so");
install("symplectic_reduction","GLD&","symplectic_reduction","./hyperellperiods.so");
install("hyperellperiods","GD0,L,p","hyperellperiods","./hyperellperiods.so");
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

/* check lattice */
/* test if Omega is symetric to current prec */
Riemann_symmetry(tau) = {
  my(Sym);
  Sym = abs(tau - mattranspose(tau));
  if(exponent(trace(Sym*mattranspose(Sym)))>-90,
      print("## Riemann_symmetry violated : Omega^t != Omega");
      return(0),
      return(1)
  );
}
/* test whether Im(Omega) > 0 */
/* this is true if and only if
   IOmega = M^t*M
   suppose Omega already tested for symmetry
 */
Riemann_positivity(tau,n=5) = {
  my(Itau=imag(tau));
  g = matsize(Itau)[1];
  for(i=1,g,if(Itau[i,i]<=0,
           printf("tau[%i,%i]<=0",i,i);
   return(0)));
  for(i=1,g,
      for(j=i+1,g,
        if(Itau[i,i]<=abs(Itau[i,j]),
           printf("tau[%i,%i]<=|tau[%i,%i]|",i,i,i,j);
           return(0));
        if(Itau[i,i]+Itau[j,j]<=2*abs(real(Itau[i,j])),
           printf("tau[%i,%i]+tau[%i,%i]<=2|Re tau[%i,%i]|",
           i,i,j,j,i,j);
           return(0));
        );
     );
  /* plus some random tests */
  X = vector(g,i,tan(3.14*random(1.)-1.55));
  for(k=1,n,
      if(X*Itau*X~<0,
        printf("X*tau*X~<0 for X=%s",X);
       return(0));
      for(i=1,g,X[i] = tan(3.14*random(1.)-1.55););
     );
  return(1)
} 


test_Riemannrelations(X) = {
  my(tau);
  /* test if the small period matrix is symetric */
  tau = hyperellperiods(X);
  if(!Riemann_symmetry(tau),error("symmetry : ",X));
  if(!Riemann_positivity(tau),error("positivity : ",X));
  return(1);
}

{
/* periods */
for(d=2,7,\\30,
  for(i=1,5,
    my(pol = random(1.*x^d));
    test_Riemannrelations(pol);
    ));
}
