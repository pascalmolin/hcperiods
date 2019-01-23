read("install.gp");
read("periods.gp");
v = readvec("genus2_100.gp");
pol=genus2red([v[1][3],v[1][2]])[3];
compare(v) = {
  if(type(v)=="t_VEC"&&#v==3,pol=[v[3],v[2]],pol=v);
  if(type(pol)=="t_VEC"&&#pol==2,pol=genus2red(pol)[3]);
  x1 = hcInit(pol);
  x2 = hcinit(pol);
  t1 = hcSmallperiods(x1);
  t2 = hyperellperiods(x2);

  /* check intersections */
  X = x2[2];
  K = x2[4][3]; \\ int
  E = [ e[2] | e <-  x2[4][1] ]; \\ edges
  F = [ y[2] | y <- x2[4][2] ]; \\ ydata
  for(i=1,3,
    for(j=i+1,4,
      [ia,ib]=E[i];
      [ic,id]=E[j];
      if(ib==ic,
        a = X[ia]; b = X[ib]; d = X[id];
        tau = arg((d-b)/(b-a));
        fb = F[i][6];
        fc = F[j][5];
        arg1 = -(Pi+tau)/2;
        arg2 =  (Pi-tau)/2;
        printf("intersection_abbd %i->%i->%i\n",ia,ib,id);
        print("  fab=",fb);
        print("  fbd=",fc);
        print("  tau=",tau);
        \\print(" roots ab ",F[i][1]);
        print([arg(fb/fc),arg1,arg2]);
      ,ia==ic,
        a = X[ia]; b = X[ib]; d = X[id];
        tau = arg((d-a)/(b-a));
        fa = F[i][5];
        fc = F[j][5];
        arg1 = -tau/2;
        arg2 = if(tau>0,Pi-tau/2,-Pi-tau/2);
        printf("intersection_abad %i->%i, %i->%i\n",ia,ib,ia,id);
        print("  fab=",fb);
        print("  fad=",fc);
        print("  tau=",tau);
        print([arg(fa/fc),arg1,arg2]);
      );
    ));

  \\ trees
  tree1 = x1[iTree];
  tree2 = apply(Vec,E);

  \\ intersections
  K1 = x1[iIntersection];
  K2 = x2[4][4];

  \\ integrals
  Int1 = 2*x1[iIntegrals][1,];
  Int2 = x2[4][3][1,];

  /* check positive */
  it1 = imag(t1);
  it2 = imag(t2);
  for(i=1,80,
    z = vector(#t1,j,random(2.)-1);
    if(z*it1*z~<0,error("positivity",z));
    if(z*it2*z~<0,error("positivity",z));
  );

  if(#t1==2,
  print("reduce X1");
  [a1,f1] = reduce_to_F2(t1);
  print("reduce X2");
  [a2,f2] = reduce_to_F2(t2);
  bestappr(f1,1000)==bestappr(f2,1000);
  );
}
