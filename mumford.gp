/* A = [a_1,a_2,.. a_{2g-1}]
   this allows to define standard loops A_i and B_i. Do not use a spanning tree
but standard loops A_i, B_i' */
hcInit_Mumford(A,provenbounds=0) = {
  my(hcStruct = vector(iMax));
  my(Dcalc,g,tree,tau,h,npoints,IntFactor,IntPoints,AB,Intersection,P);
  my(Omega0,Omega1);
  g = floor((#A-1)/2);
  hcStruct[iRoots] = A;
  hcStruct[iGenus] = g;  
  hcStruct[iIntDcalc] = Dcalc = precision(1.);
  /* low precision parameters */
  default(realprecision,19);
  gettime();
  /* this changes */
  my(tree=vector(2*g));
   for(i=1,g,tree[i] = [2*i-1,2*i];tree[g+i] = [2*i,2*i+1]);
   for(i=1,2*g,if(real(A[tree[i][1]])>real(A[tree[i][2]]),t = tree[i][1];
tree[i][1]=tree[i][2]; tree[i][2]=t));
  hcStruct[iTree] =  tree;
  hcStruct[iIntTau] = tau = Pi/8; /* false but does not matter */
  msgtimer("tree");
  /* prepare integration */
  my(D = Dcalc*log(10));
  tmp = integration_parameters(A,tree,tau,D,provenbounds);
  hcStruct[iIntH] = h = precision(tmp[1],Dcalc);
  hcStruct[iIntNpoints] = npoints = tmp[2];
  msgtimer("bounds");
  /* return to normal precision */
  default(realprecision,Dcalc);
  tmp = integration_points_thsh(h,npoints); \\ bad...
  IntFactor = tmp[1];
  IntPoints = tmp[2];
  hcStruct[iIntFactor] = IntFactor;
  hcStruct[iIntPoints] = IntPoints;
  msgtimer("integration points");
  /* compute big period matrix of periods C_1,..C_2g */
  coh1x_homC = periods_spanning(hcStruct);
  msgtimer("periods");
  /* intersection matrix of basis C_1,.. C_2g */
  IntC = intersection_spanning(A,tree);
  hcStruct[iIntersection] = IntC;
  msgtimer("intersection");
  /* column base change from (A_i,B_i) to C_i */
  /* that is, P^t * Intersection * P = J_g */
  ABtoC = matid(2*g);
  \\ABtoC = symplectic_reduction(IntC);
  hcStruct[iABtoC] = ABtoC;
  /* change basis */
  /* put everything in the struct */
  coh1x_homAB = coh1x_homC * ABtoC; /* now matrix in (A,B) basis */
  Omega0 = vecextract(coh1x_homAB,Str("1..",g),Str("1..",g));
  Omega1 = vecextract(coh1x_homAB,Str("1..",g),Str(g+1,".."));
  hcStruct[iOmega0] = Omega0; \\ on A basis
  hcStruct[iOmega1] = Omega1; \\ on B basis
  /* with dual basis for H^1 */
  hcStruct[iTau] = tau = Omega1^(-1)*Omega0; \\ van Wamelen convention
  /* init Abel-Jacobi map */
  \\hcStruct[iP0] = P0 = tree[1][1]; /* choose central root */
  hcStruct[iP0] = P0 = 1; /* choose first root */
  \\hcStruct[iAJRoots] = (Omega1)^(-1)*coh1x_homAB*AJWeierstrass(hcStruct);
  hcStruct[iAJRoots] = AJWeierstrass(hcStruct);
  /* reduction map to (A,B) lattice base */
  /* choose basis [Id,tau] */
  \\RLleft = concat(matid(g),matrix(g,g))~;
  \\RLright = concat(real(tau)~,imag(tau)~)~;
  \\hcStruct[iReduce] = (concat(RLleft,RLright))^(-1);
  /* choose Omega0,Omega1 */
  M = matrix(2*g,2*g);
  for(i=1,g,for(j=1,g,
    M[  i,  j] = real(Omega0[i,j]);
    M[g+i,  j] = imag(Omega0[i,j]);
    M[  i,g+j] = real(Omega1[i,j]);
    M[g+i,g+j] = imag(Omega1[i,j]);
  ));
  hcStruct[iReduce] = M^(-1);
  /* it's done */
  return(hcStruct);
}

