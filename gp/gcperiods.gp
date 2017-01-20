/**

======================================
Abel-Jacobi map on hyperelliptic curves
=======================================

computing period matrix of hyperelliptic curves

*/

msgtimer(s,level=1) = {
  if(default(debug)>=level,printf("  Time %s : %ld\n",s,gettime()));
}
msgdebug(s,header="  ",level=1) = {
  if(default(debug)>=level,printf("%s%s\n",header,s));
}
/**

Numerical integration with Gauss Chebychev
------------------------------------------
*/
/* compute prod sqrt(z-u_i) as sqrt prod(z-u_i)
   track the number of loops around 0
   while multiplying the (z-ui).
 */
/* sign of the imaginary part with positive zero */
isign(x) = { my(s);s=sign(imag(x));if(s,s,1);}
sqrt_affinereduction_gc(u,z,flag=0) = {
  my(z1,i1,z2,i2,sgn);
  \\ if we want the long method, or there is only one root
  flag=1; /* FIXME! */
  if(flag||#u==1,
    prod(i=1,#u, if(real(u[i])>0,I*sqrt(u[i]-z),sqrt(z-u[i])))
    , \\ otherwise
    sgn = 1; \\ track the number of loops
    z1 = z-u[1];
    /* if the point is on the right, take non-standard cut */
    if(imag(z1)==0&&real(z1)<0,z1=-z1;sgn*=-I);
    i1 = isign(z1);
    for(i=2,#u,
      z2 = z-u[i];
      if(imag(z2)==0&&real(z2)<0,z2=-z2;sgn *=-I);
      i2 = isign(z2);
      z1 *= z2;
      if(i1==i2,
        i1 = isign(z1);
        if(i1!=i2,sgn*=-1),
        i1 = isign(z1)
        ));
    sqrt(z1)*sgn
    );
}
/*
 compute the column of period integrals
 [ int x^k/y dx ]
*/
integrals_gc(u, g, n, flag=0) = {
  my(res = vector(g));
  for(k=1,n,
    my(xk = cos((2*k-1)*Pi/(2*n)));
    my(t = 1/sqrt_affinereduction_gc(u,xk,flag));
    res[1] += t;
    for(i=2, g, t *= xk; res[i] += t);
    );
  res * Pi/n;
}
/* make u vectors : remove two roots */
ab_points_gc(A, ia, ib) = {
    my(a = A[ia]); my(b=A[ib]); my(j=0);
    my(u=vector(#A-2));
    for(i=1,#A,i==ia||i==ib||u[j++]=(2*A[i]-a-b)/(b-a));
    u;
}
polshift(v,c)=for(i=2,#v,forstep(j=#v,i,-1,v[j]+=c*v[j-1]));v;
integrals_edge_gc(A, g, i1, i2, flag=0) = {
  my(u = ab_points_gc(A, i1, i2));
  msgdebug([A[i1],A[i2]], "=========== gc ==========\n    a,b =", 2);
  my(n = gc_n(u));
  my(res = integrals_gc(u, g, n));
  msgdebug(res, "    int =", 2);
  /* shift */
  my(a,b,c);
  a = A[i1]; b = A[i2];
  c = (b+a)/(b-a);
  res = polshift(res, c);
  msgdebug(res, "    shifted =", 2);
  /* fix factors */
  my(geom_factor = (b-a)/2);
  my(fact = I/sqrt(geom_factor)^#u);
  res[1] *= fact;
  for(j=2,g,
          fact *= geom_factor;
          res[j] *= fact;
     );
  msgdebug(res, "    factors =", 3);
  \\ to track signs, return also the values of the sqrt chosen at end points
  vplus = sqrt_affinereduction_gc(u,1,flag)*sqrt(geom_factor)^#u/geom_factor;
  vmoins = sqrt_affinereduction_gc(u,-1,flag)*sqrt(geom_factor)^#u/geom_factor;
  return([res,vmoins,vplus]);
}
/** == compute integration parameters == */
dist_ellipse(p,a) = {
    my(b,t);
    my(x,y);
    my(eps=0.*1e8);
    b = sqrt(a^2-1);
    x=abs(real(p));y=abs(imag(p));
    if(x==0,abs(b-y),
            y==0,abs(a-x),
            t = atan(a*y/(b*x));
            ft = 1;
            my(i);
            while(abs(ft) > eps,
                \\if(i++ % 8 == 0,print(i,[p, a ,t]));
                st = sin(t);
                ct = cos(t);
                ft = ct*st-x*a*st+y*b*ct;
                fpt = ct^2-st^2-x*a*ct-y*b*st;
                t -= ft / fpt;
                );
            abs(p-a*cos(t)-I*b*sin(t))
      );
}
gc_n(u,prec,mult=.25) = {
    \\my(A,B);
    if(!prec,prec=precision(1.));
    au = [ abs(z-1)+abs(z+1) | z<-u ] / 2;
    v = vecsort(au, , 1);
    a = au[v[1]];
    m = 1; while(m<#u&&au[v[m+1]]-a<.001,m++);
    r = a + sqrt(a^2-1);
    /* first estimate */
    M = 1/sqrt(a^m*prod(i=m+1,#u,dist_ellipse(u[v[i]],a)));
    A = prec*log(10) + log(2*Pi*M) + 1;
    B = 2*log(r);
    msgdebug(A/B,"n1 = ",2);
    n = ceil(A / B);
    /* second */
    a = a*(1-mult/n);
    M = 1/sqrt(prod(i=1,#u,dist_ellipse(u[i],a)));
    A = prec*log(10) + log(2*Pi*M) + 1;
    r = a + sqrt(a^2-1);
    B = 2*log(r);
    msgdebug(A/B,"n2 = ",2);
    ceil(A/B);
}
r_3(a,b,c) = {
  (abs(c-a)+abs(c-b))/abs(b-a);
}
r_edge(A,i,j) = {
  rij = oo;
  for(k=1,#A,
    if(k==i||k==j,next());
      rij = min(rij,r_3(A[i],A[j],A[k]));
      );
   return(rij);
}
r_tree(A,edges) = {
  my(r = oo);
  for(k=1,#edges,
    r = min(r,r_edge(A,edges[k][1],edges[k][2]));
  );
  return(r);
}
/* *************************************************************************
 *  period matrix
 * ************************************************************************* */

/* *************************************************************************
 *  step 1: maximal-flow spanning tree
 * ************************************************************************* */
max_spanning_gc(A,nedges=0)= {
  /* low precision parameters */
  my(n,r_v,z,taken,tree,rmin);
  localprec(19);
  n = #A;
  if(n<=0,error(A));
  if(!nedges,nedges=n-1); /* spanning tree */
  r_v=vector(n*(n-1)/2);
  /* build list of edges ordered by r */
  z=0;
  for(i=1,n,
    for(j=i+1,n,
      r_v[z++] = [r_edge(A,i,j),i,j];
      );
    );
  r_v = vecsort(r_v,1); \\ sort according to r
  taken = vector(n,i,0); \\ mark vertex if in tree
  tree = vector(nedges);
  /* this is not optimal... too many loops */
  z=n*(n-1)/2; \\select edge with maximum r
  rmin = r_v[z][1];
  for(k=1,nedges,
    z=n*(n-1)/2;
    /* consider next edge with exactly one vertex taken (no cycle) */
    while(k>1 && taken[r_v[z][2]]==taken[r_v[z][3]],z--);
    /* this is the best edge allowed */
    i = r_v[z][2]; j = r_v[z][3];
    /* ensure i already in tree */
    if(taken[j],[i,j]=[j,i]);
    tree[k] = [i,j];
    rmin = min(rmin,r_v[z][1]);
    taken[i] = taken[j] = 1;
    );
  return([tree,rmin]);
}

/* *************************************************************************
 *  step 2 : intersection matrix of selected edges
 * ************************************************************************* */

/* end intersection I[ab].I[bd] */
intersection_abbd_gc(A,ia,ib,ic,id) = {
    my(k,a,b,d,tau,Aprim,fad,j);
    k = #A;
    a = A[ia]; b = A[ib]; d = A[id];
    tau = arg((d-b)/(b-a)); 
    /* \frac{\sqrt{d-b}^{k}\prod_{λ\neq b,d}\sqrt{x_β(b)-ψ_β(λ)}} */
    \\fbd = sqrt((d-b))^k
        \\*prod(i=1,k,if(i==ib||i==id,1,sqrt(-1-(2*A[i]-b-d)/(d-b))));
    u = ab_points_gc(A,ib,id);
    fbd = sqrt((d-b))^k*sqrt_affinereduction_gc(u,-1);
    \\fab = sqrt((b-a))^k
        \\*prod(i=1,k,if(i==ib||i==ia,1,sqrt(+1-(2*A[i]-b-a)/(b-a))));
    u = ab_points_gc(A,ia,ib);
    fab = sqrt((b-a))^k*sqrt_affinereduction_gc(u,1);
    if(default(debug)>=1,
            endarg = arg(fab/fbd);
            arg1 = -(Pi+tau)/2;
            arg2 = (Pi-tau)/2;
            if(abs(endarg-arg1)>0.0001&&abs(endarg-arg2)>0.0001,
                print("intersection_abbd");
                print("  tau=",tau);
                printf("  arg[fab/fbd]=%1.4f should be %1.4f or %1.4f\n",endarg,arg1, arg2);
                printf("tau/2=%1.3f ; 2*tau=%1.3f ; (Pi-tau)/2=%1.3f ; Pi-2*tau=%1.3f ; \n",
                    tau/2,2*tau,(Pi-tau)/2,Pi-2*tau);
                error();
              );
      );
    if(imag(fbd/fab)<0,return(-1),return(1));
}
intersection_abad_gc(A,ia,ib,ic,id) = {
    my(k,a,b,d,tau,u,fad,j);
    k = #A;
    a = A[ia]; b = A[ib]; d = A[id];
    tau = arg((d-a)/(b-a)); 
    \\ \frac{\sqrt{d-b}^{k}\prod_{λ\neq b,d}\sqrt{x_β(b)-ψ_β(λ)}}
    \\fad = sqrt((d-a))^k
        \\*prod(i=1,k,if(i==ia||i==id,1,sqrt(-1-(2*A[i]-a-d)/(d-a))));
    u = ab_points_gc(A,ia,id);
    fad = sqrt((d-a))^k*sqrt_affinereduction_gc(u,-1);
    \\fab = sqrt((b-a))^k
        \\*prod(i=1,k,if(i==ib||i==ia,1,sqrt(-1-(2*A[i]-b-a)/(b-a))));
    u = ab_points_gc(A,ia,ib);
    fab = sqrt((b-a))^k*sqrt_affinereduction_gc(u,-1);
    if(default(debug)>=1,
            endarg = arg(fab/fad);
            arg1 = -tau/2;
            arg2 = if(tau>0,Pi-tau/2,-Pi-tau/2);
            if(abs(endarg-arg1)>0.0001&&abs(endarg-arg2)>0.0001,
                print("intersection_abad");
                print("  tau=",tau);
                printf("  arg[fab/fad]=%1.4f should be %1.4f or %1.4f\n",endarg,arg1, arg2);
                printf("tau/2=%1.3f ; 2*tau=%1.3f ; (Pi-tau)/2=%1.3f ; Pi-2*tau=%1.3f ; Pi-tau/2=%1.3f ; Pi+tau/2=%1.3f \n", tau/2,2*tau,(Pi-tau)/2,Pi-2*tau,Pi-tau/2,Pi+tau/2);
                error();
              );
      );
    if(imag(fab/fad)<0,return(-1),return(1));
}
intersection_abcb_gc(A,ia,ib,ic,id) = {
    my(k,a,b,c,tau,Aprim,fcb,fab,j);
    k = #A;
    a = A[ia]; b = A[ib]; c = A[ic];
    tau = arg((b-c)/(b-a)); 
    \\ \frac{\sqrt{d-b}^{k}\prod_{λ\neq b,d}\sqrt{x_β(b)-ψ_β(λ)}}
    \\fcb = sqrt((b-c))^k
        \\*prod(i=1,k,if(i==ic||i==ib,1,sqrt(+1-(2*A[i]-c-b)/(b-c))));
    Aprim = ab_points_gc(A,ic,ib);
    fcb = sqrt((b-c))^k * sqrt_affinereduction_gc(Aprim,+1);
    \\fab = sqrt((b-a))^k
        \\*prod(i=1,k,if(i==ib||i==ia,1,sqrt(+1-(2*A[i]-b-a)/(b-a))));
    Aprim = ab_points_gc(A,ia,ib);
    fab = sqrt((b-a))^k*sqrt_affinereduction_gc(Aprim,+1);
    if(default(debug)>=1,
            endarg = arg(fab/fcb);
            \\arg1 = tau/2;
            arg1 = if(tau>0,Pi-tau/2,-Pi-tau/2);
            \\arg2 = Pi-tau/2;
            arg2 = -tau/2;
            if(abs(endarg-arg1)>0.0001&&abs(endarg-arg2)>0.0001,
                print("intersection_abcb");
                print("  tau=",tau);
                printf("  arg[fab/fcb]=%1.4f should be %1.4f or %1.4f\n",endarg,arg1, arg2);
                printf("tau/2=%1.3f ; 2*tau=%1.3f ; (Pi-tau)/2=%1.3f ; Pi-2*tau=%1.3f ; Pi-tau/2=%1.3f ; Pi+tau/2=%1.3f \n", tau/2,2*tau,(Pi-tau)/2,Pi-2*tau,Pi-tau/2,Pi+tau/2);
                error();
              );
      );
    if(imag(fab/fcb)<0,return(-1),return(1));
}
/* TODO: check aligned points ... */
intersection_gc(A,ia,ib,ic,id) = {
    \\ do we have a common point ?
        if(ia==ib||ic==id, 0,
                (ia==ic&&ib==id)||(ia==id&&ib==ic), 0,\\ self intersection
                ia!=ic&&ia!=id&&ib!=ic&&ib!=id,0,\\ no intersection
                ia==ic, +intersection_abad_gc(A,ia,ib,ic,id),
                ib==ic, +intersection_abbd_gc(A,ia,ib,ic,id),
                ia==id, -intersection_abbd_gc(A,ic,id,ia,ib),
                ib==id, +intersection_abcb_gc(A,ia,ib,ic,id),
                error("this intersection should not occur",[ia,ib,ic,id]);
          );
}
intersection_spanning_gc(A,tree) = {
    my(n,res,s);
    n = #tree;
    res = matrix(n,n);
    for(i=1,#tree,
            for(j=i+1,#tree,
                s = intersection_gc(A,tree[i][1],tree[i][2],tree[j][1],tree[j][2]);
                res[i,j] =  s;
                res[j,i] = -s;
               );
       );
    res;
}
/* *************************************************************************
 *  step 3 : period matrix of edges + intersection matrix
 * ************************************************************************* */
periods_spanning_gc(A,tree,g) = {
  if(1, /* forget about end points */
  Mat(vector(#tree,k,integrals_edge_gc(A,g,tree[k][1],tree[k][2])[1]~));
  ,
  res = vector(#tree);
  endpointsvalues=vector(#tree);
  for(k=1,#tree,
    tmp = integrals_edge_gc(A,g,tree[k][1],tree[k][2]);
    /* should use final values to find intersection if sure that no
     * crossing edges TODO */
    res[k] = tmp[1]~;
    endpointsvalues[k] = [tmp[2],tmp[3]];
  );
  Mat(res);
  );
}

/* *************************************************************************
 *  step 4 : retrieve symplectic basis
 * ************************************************************************* */

/* assume M is invertible, antisymetric and of even dimension 2g */
symplectic_reduction(M) = {
  my(g=matsize(M)[1]);
  if(g%2,return(0),g/=2);
  my(N = Vec(M)); \\ columns matrices
  my(P = Vec(matid(2*g)));
  my(i,j,ii,jj);
  /* algorithm similar to snf
     take minimal N[i,j] non-zero value above diagonal
     1. transform it to be the gcd over the line
     2. reduce lines i and j to have symplectic
        [0,1;-1,0] and zeroes outside
     3. move these in front
  the same operations are done on P to track
  the base change.
  */
  \\ find a non zero value;
  \\ TODO better to take minimum
    /* main loop on symplectic 2-subspaces */
    for(d = 1, g,
        if(1,
          my(res = Mat(vector(2*g,k,if(k<=g,P[2*k-1],P[2*(k-g)]))));
          \\print(res)
          \\print(Mat(N)); \\ FIXME
          );
  
        /* find a pivot value */
        i = 2*d-1; j = i+1;
        while(!N[j][i],
            print("i",i,"j",j,"->",N[j][i]);
            if(j++>2*g,j=i+++1);
          );
        \\ i and j give a and b
        if(N[j][i]<0,
           P[j] = -P[j]; 
           N[j] = -N[j]; 
           for(l=1,2*g,N[l][j] = -N[l][j]); 
        );
        ii=i; jj=j;
        while(jj,
        \\print("d=",d," -> i=",i,"&j=",j);
        /* remove all coeffs on column j */
          for(k=2*d-1,2*g, \\ the column
            \\print("column ",k);
            if(k==jj||!N[k][ii],next()); \\ next if already zero
            \\print(Mat(N));
            if(N[k][ii]%N[jj][ii]==0,
              /* N[jj][ii] divides N[k][ii], set it to zero */
              my(q = N[k][ii]/N[jj][ii]);
              P[k] = P[k]-q*P[jj];
              N[k] = N[k]-q*N[jj];
              for(l=1,2*g,N[l][k] = N[l][k] - q*N[l][jj];);
              ,
              /* otherwise replace N[jj][ii] by the gcd */
              tmp = bezout(N[jj][ii],N[jj][k]);
              u = tmp[1]; v = tmp[2]; d = tmp[3];
              P[ii] = u*P[ii]+v*P[jj];
              N[ii] = u*N[ii]+v*N[jj];
              for(l=1,2*g,N[l][ii] = u*N[l][ii] + v*N[l][jj];);
              \\ now resume on N[jj][ii]
              k--
              );
            );
          if(jj>ii,ii=j;jj=i,jj=0);
  );
  /* and put the 2-2 symplectic at 2*d-1, 2*d */
  tmpN = N[2*d-1]; tmpP = P[2*d-1];
  N[2*d-1] = N[i]; P[2*d-1] = P[i];
  N[i] = tmpN; P[i] = tmpP;
  tmpN = N[2*d]; tmpP = P[2*d];
  N[2*d] = N[j]; P[2*d] = P[j];
  N[j] = tmpN; P[j]= tmpP;
  );
  /* and finally canonical basis */
  res = Mat(vector(2*g,k,if(k<=g,P[2*k-1],P[2*(k-g)])));
  \\print(mattranspose(res)*M*res);
  return(res); 
}
/* inverse of a symplectic matrix */
symplectic_inverse(M) = {
  my(g = matsize(M)[1]/2);
  /* M=[A,B;C,D] -> M^(-1)=[D^t,-B^t,-C^t,A^t] */
  my(res=matrix(2*g,2*g));
  for(i=1,g,
    for(j=1,g,
      /* A block */
      res[i,j] = M[g+j,g+i];
      /* D block */
      res[g+i,g+j] = M[j,i];
      /* B block */
      res[i,g+j] = -M[g+j,i];
      /* C block */
      res[g+i,j] = -M[j,g+i];
   );
 );
return(res);
}

/* *************************************************************************
 *  put everything together
 * ************************************************************************* */

/** Remark on matrices conventions :
 during computations, the tree corresponds to a homology basis C_1,..C_2g ; we
have the canonical symplectic homology basis A_1,..A_g,B_1,..B_g. We compute
first the period matrix C whose columns are the periods of the canonical
differential forms along the paths C_j :

coh1x_homC =    C_1,   C_2g
               +-----------+
             1 |
             x |
 
Then we
compute the intersection matrix of C_j

IntC = C_1,   C_2g
     +-----------+
 C_1 |
 C_2 |

and we perform symplectic reduction to have the base change to a canonical
basis (which is not Mumford's in general)

ABtoC = A_1,  A_g,B_1  B_g
     +------------------+
C_1  |
...  |
C_2g |

satisfying :

(ABtoC)^t * IntC * ABtoC = J_g = intersection in (A,B) basis.

Now we have the big period matrix of canonical differentials in symplectic basis

coh1x_homAB = coh1x_homC*ABtoC = [ Omega_0 Omega_1 ]
    =      A_1,.. B_g
          +----------+
       1  |
      ... |
     x^g-1|

and we take the (A)-dual cohomology basis to have a Siegel tau matrix

cohWto1x = 1,x,.. x^g-1
          +------------+
       w_1|
       ...|
       w_g|

cohWto1x*Omega_0 = Id
cohWto1x*Omega_1 = tau = Omega_0^(-1)*Omega_1

As for Abel-Jacobi map : we compute it in the 1,x cohomology basis
and return its left-product by cohWto1x.

For Weierstrass points Pi, we compute half loops P0Pi in C_i basis

P0Pi = P0P1,... P0P2g+1  (one column is zero for Pi=P0)
      +--------------+
  C_1 |
  ... |

then in A,B basis

(ABtoC)^t * P0Pi

and finally return the integration

cohWto1x * coh1x_homAB * (ABtoC)^t * P0Pi
*/

/**
the period matrix is computed in some symplectic homology basis
determined during computations.
*/
hcperiods(f,flag) = {
  my(A);
  if(type(f)=="t_VEC",A=f,
     type(f)=="t_COL",A=Vec(f),
     type(f)=="t_POL",A=polroots(f),
     error("unrecognized type for f ", f);
     );
  my(g,tree,tmp,ABtoC,IntC);
  g = floor((#A-1)/2);
  gettime();
  tmp = max_spanning_gc(A,2*g);
  tmp = max_spanning(A,2*g);
  msgtimer("tree");
  tree = tmp[1];
  msgdebug(tree,"  ",2);
  /* compute big period matrix of periods C_1,..C_2g */
  my(coh1x_homC = periods_spanning_gc(A,tree,g));
  msgtimer("periods");
  msgdebug(coh1x_homC,"  ",2);
  /* intersection matrix of basis C_1,.. C_2g */
  IntC = intersection_spanning_gc(A,tree);
  msgtimer("intersection");
  msgdebug(IntC,"  ",2);
  /* column base change from (A_i,B_i) to C_i */
  /* that is, P^t * Intersection * P = J_g */
  ABtoC = symplectic_reduction(IntC);
  /* change basis */
  /* put everything in the struct */
  coh1x_homAB = coh1x_homC * ABtoC; /* now matrix in (A,B) basis */
  my(Omega0,Omega1);
  Omega0 = vecextract(coh1x_homAB,Str("1..",g),Str("1..",g));
  Omega1 = vecextract(coh1x_homAB,Str("1..",g),Str(g+1,".."));
  if(flag==0, Omega1^(-1)*Omega0, /* default, tau matrix */
     flag==1, concat(Omega0,Omega1),
     error("unknown flag"));
}

/* *************************************************************************
 * tests
 * ************************************************************************* */

/* Gauss's AGM
   I = int_a^b dx/sqrt((x-a)(x-b)(x-c))
     = Pi/agm(sqrt(c-a),sqrt(c-b))
*/
gauss_agm(P) = {
  my(a,b,c);
  [a,b,c] = P;
  /* a>b>c */
  if(0,
  Mat([Pi/agm(a,b),I*Pi/agm(a,c)]);
  ,
  Mat([I*Pi/agm(sqrt(b-c),sqrt(a-c)),Pi/agm(sqrt(a-b),sqrt(a-c))]);
  );
}

/* Richelot's algorithm
   Iu=int_u^u' S(x)/sqrt(-P(x)) dx
   u<u'<v<v'<w<w'
   P(x)=(x-u)(x-u')(x-v)(x-v')(x-w)(x-w')
   S(x)=lx+m
   P=[u,u',v,v',w,w']
   S=[l,m]
*/
richelot(P,S)= {
  my(a,b,c,ap,bp,cp,D,A,B,C,dba,dcb,dca,a1,a1p,b1,b1p,c1,c1p,D1);
  a=P[1];ap=P[2];b=P[3];bp=P[4];c=P[5];cp=P[6];
  D=1;
  prec = default(realprecision);
  localprec(prec+5);
  epsilon=10^-(default(realprecision));
  while(1,
    print([a,ap,b,bp,c,cp]);
    A=sqrt((c-b)*(c-bp)*(cp-b)*(cp-bp));
    B=sqrt((c-a)*(c-ap)*(cp-a)*(cp-ap));
    C=sqrt((b-a)*(b-ap)*(bp-a)*(bp-ap));
    dba = b+bp-a-ap;
    dcb = c+cp-b-bp;
    dca = c+cp-a-ap;
    a *= ap; b *= bp; c *= cp;
    a1 = (c-a-B)/dca;
    a1p= (b-a-C)/dba;
    b1 = (b-a+C)/dba;
    b1p= (c-b-A)/dcb;
    c1 = (c-b+A)/dcb;
    c1p= (c-a+B)/dca;
    D1 = 4*D*(a*dcb-b*dca+c*dba)/(dba*dcb*dca);
    a=a1; ap=a1p; b=b1; bp=b1p; c=c1; cp=c1p; D=D1;
    /* test */
    if(vecmax(abs([a,b,c]-[ap,bp,cp]))<epsilon,
			return(Pi*sqrt(D)*(S[1]*a+S[2])/((a-b)*(a-c)))
    )
  );
}
/* 0.4357814509980649277362248425469236 */
intnum(x=[1,-.5],[2,-.5],1/sqrt(-prod(j=1,6,x-j)))
richelot([1..6],[0,1])
/* 0.6969159933749301638525527624477020 */
intnum(x=[1,-.5],[2,-.5],x/sqrt(-prod(j=1,6,x-j)))
richelot([1..6],[1,0])

X1 = [-1,I,-I,1,1+I]; \\ house
intnum(x=[-1,-.5],[1,-.5],1/sqrt(-prod(j=1,#X1,x-X1[j])))

X2 = [-1,0,-I,1-I,I,1+I];
intnum(x=[-1,-.5],[0,-.5],1/sqrt(-prod(j=1,#X2,x-X2[j])))







\\step(z) = if(z,z,1);
\\checkint() = {
\\  A = polroots(x^7-x-1);
\\  u = ab_points_gc(A,6,7);
\\  a = A[6]; b = A[7];
\\  ref = intnum(x=[-1,-.5],[1,-.5],[1,x,x^2]/(sqrt(1-x^2)*sqrt_affinereduction_gc(u,x)));
\\  r1 = integrals_gc(u, 3, 100);
\\  r2 = sqrt(I)*intnum(x=[a,-.5],[b,-.5],[1,x,x^2]/step(sqrt_affinereduction_gc(A,x)));
\\  print(ref);
\\  print(r1);
\\  print(r2);
\\}


/* check if same as older DE code */
check(f) = {
  my(c,tree,g,p1,p2);
  c = hcInit(f);
  tree = c[iTree];
  g = c[iGenus];
  p1 = Vec(periods_spanning(c));
  p2 = Vec(periods_spanning_gc(c[iRoots],tree,g));
  for (i=1,#p1,
    if(vecsum(abs(p1[i]-p2[i]))>.0000001,
    my(u = ab_points_gc(c[iRoots], tree[i][1], tree[i][2]));
    print("i=",i);
    print("a,b=",[c[iRoots][tree[i][1]],c[iRoots][tree[i][2]]]);
    print("u=",u);
    print(p1[i]);
    print(p2[i]);
    ));
}
/* test if Omega is symetric to current prec */
Riemann_symmetry(Omega) = {
  my(Sym);
  Sym = abs(Omega - mattranspose(Omega));
  if(trace(Sym*mattranspose(Sym))>1e-20,
      print("## Riemann_symmetry violated : Omega^t != Omega");0,1);
}
/* test whether Im(Omega) > 0 */
/* this is true if and only if
   IOmega = M^t*M
   suppose Omega already tested for symmetry
 */
Riemann_positivity(Omega,n=5) = {
  my(IOmega,g,X);
  IOmega=imag(Omega);
  g = matsize(IOmega)[1];
  for(i=1,g,if(IOmega[i,i]<=0,
           msgdebug(Strprintf("tau[%i,%i]<=0",i,i),"  *** ");
   return(0)));
  for(i=1,g,
      for(j=i+1,g,
        if(IOmega[i,i]<=abs(IOmega[i,j]),
           msgdebug(Strprintf("tau[%i,%i]<=|tau[%i,%i]|",i,i,i,j),"  *** ");
           return(0));
        if(IOmega[i,i]+IOmega[j,j]<=2*abs(real(IOmega[i,j])),
           msgdebug(Strprintf("tau[%i,%i]+tau[%i,%i]<=2|Re tau[%i,%i]|",
           i,i,j,j,i,j),"  *** ");
           return(0));
        );
     );
  /* plus some random tests */
  X = vector(g,i,tan(3.14*random(1.)-1.55));
  for(k=1,n,
      if(X*IOmega*X~<0,
        msgdebug(Str("X*tau*X~<0 for X=",X),"  *** ");
       return(0));
      for(i=1,g,X[i] = tan(3.14*random(1.)-1.55););
     );
  return(1)
} 
test_Riemannrelations(pts) = {
  my(Omega = hcperiods(pts));
  /* test if the small period matrix is symetric */
  if(!Riemann_symmetry(Omega),error("symmetry : ",pts));
  if(!Riemann_positivity(Omega),error("positivity : ",pts));
  return(1);
}
ptsRandom(g=2) = {
  vector(2*g+1,k,tan(3.1*(random(1.)-.5))+I*tan(3.1*(random(1.)-.5)));
}
test_Riemannrelations_random(n=100) = {
  my(pts);
  for(i=1,n,
    pts = ptsRandom(random(7)+1);
    test_Riemannrelations(pts) || return(0);
    );
  return(1);
}

/* period([a,b,c]) ~ agm()
   period([a,b,c,d,e,f])~richelot([a,b,c,d,e,f],[0,1])
*/
gccompare(P) = {
  my(g,p,i);
  if(type(P)=="t_POL",P=polroots(P));
  if(#P==3,
    gettime();
    g = gauss_agm(P);
    print("temps agm : ",gettime(),"ms.");
    g = gauss_agm(P);
    print(g);
    print("temps agm : ",gettime(),"ms.");
    p = hcperiods(P,1);
    print("temps integration : ",gettime(),"ms.");
    print(p);
    print("erreur : ",g-p);
    g,
    #P==6,
    gettime();
    g = richelot(P,[0,1]);
    print("temps richelot : ",gettime(),"ms.");
    p = hcperiods(P,1);
    print("temps integration : ",gettime(),"ms.");
    print("erreur : ",g-I*p);
    i = period_intnum(P);
    print("temps pari intnum : ",gettime(),"ms.");
    print("erreur : ",g-I*i);
    g,
    period_old(P));
}
