/**

======================================
Abel-Jacobi map on hyperelliptic curves
=======================================

computing period matrix and Abel-Jacobi map for hyperelliptic curves
using double-exponential integration between Weierstrass points.

The image of Abel-Jacobi is given in (\R/\Z)^2g.


*/

msgtimer(s,level=1) = {
  if(default(debug)>=level,printf("  Time %s : %ld\n",s,gettime()))
}
msgdebug(s,header="  ",level=1) = {
  if(default(debug)>=level,printf("%s%s\n",header,s))
}

/**

Structure
---------

We compute the following data (see differences with Magma below)

  C = hcInit(pol) : initialize the hyperelliptic curve structure : compute
    invariants for the AbelJacobi map. Returns the structure.

gp stored values :
  hcPol(C) : polynomial
  hcGenus(C) : the genus
  hcTree(C) : list of directed edges giving a basis of homology C1,.. C_{2g}
  hcIntersection(C) : matrix of intersection between edges of hcTree
  hcABtoC(C) : symplectic base change to get a symplectic basis
    from the above.
  hcSymplecticBasis(C) : a canonical homology basis
  hcBigperiods(C) : big period matrix in the canonical homology basis
  hcSmallperiods(C) : small period matrix
  hcAJRoots : AbelJacobi map on Weierstrass points
  hcP0(C) : base point for Abel-Jacobi

Computations :
  hcAbelJacobi(C,p) : image of P
  hcAbelJacobiVec(C,v) : image of a vector of points (more efficient)

internal use only :
  hcIntegrationTau : max value tau for integration
  hcIntegrationPoints

Magma stored values :

  Genus: Genus, g, of the curve for which this is the Jacobian.
  SmallPeriodMatrix: The gxg matrix in Siegel-upperhalf space
  BigPeriodMatrix: The gx2g period matrix
  HyperellipticPolynomial: Polynomial f such that y = f(x) defines the
    hyperelliptic curve
  OddDegree: Boolean value. True if f has odd degree, false otherwise
  EvenLeadingCoefficient: If f was given with even degree, it's leading
    coefficient. Otherwise undefined
  OddLeadingCoefficient: If f was given with odd degree, it's leading
    coefficient. Otherwise the leading coefficient of the polynomial obtained
    by sending InfiniteRoot to infinity.
  InfiniteRoot: If f had even degree, one of its roots. This root is sent to
    infinity in order to get a odd degree model.
  Roots: Roots of the odd degree polynomial. Either the polynomial
    given or after it was put in odd degree form. The order of these roots
    are fixed (using CompareRoots).
  EndPoints: Points in the complex plane used as endpoints of straight line
    segments used to build up a homology basis.
  AHomologyBasis: A list of 2g lists of indices into EndPoints giving closed
    loops on the Riemann surface forming a homology basis.
  ABtoC: The 2gx2g matrix with which AHomologyBasis must be
    multiplied to form a symplectic homology basis (which has a nice
    SmallPeriodMatrix).
  BasePoint: Base point for the paths. Also used as the base point from
    which the integral to infinity is computed.
  DirectedEdges: A list of pairs of indices into EndPoints. Each pair [i,j]
    represents the straight line from EndPoints[i] to EndPoints[j].
  PathIntegrals: For each entry in DirectedEdges this contains a corresponding
    entry giving the integral with respect to all g differentials along this
    straight line.
  PathSigns: Whether analytic continuation of y along the path changes the
    sign of (magma's value of) Sqrt[f(x)].
  PathGraph: The graph formed by the DirectedEdges.
  InfiniteIntegrals: The integrals with respect to the g differentials from
    BasePoint to infinity.
  ThetaNullCharacteristics: List of theta characteristics at which the
    thetanuls (at SmallPeriodMatrix) are known. Not computed here, but in
    jacmaps.m
  ThetaNullValues: The thetanull values at these characteristics
  DifferentialChangeMatrix: Change of variable matrix for the difference in
    using x^i dx/y for the even and odd degree models.

*/

iGenus            =  1 ;
iRoots            =  2 ;
iTree             =  3 ;
iIntersection     =  4 ;
iABtoC            =  5 ;
iSymplecticBasis  =  6 ;
iOmega0           =  7 ; \\ should perhaps store its inverse...
iOmega1           =  8 ;
iIntersection     =  9 ;
iIntNpoints       =  10;
iIntH             =  11;
iIntTau           =  12;
iIntDcalc         =  13;
iIntFactor        =  14;
iIntPoints        =  15;
iP0               =  16;
iAJRoots          =  17;
iReduce           =  18;
iTau              =  19;
iIntegrals        =  20;
iMax              =  iIntegrals;

/**

Numerical integration with the double-exponential formula
---------------------------------------------------------

*/
/** == compute integration parameters == */

/* limit value for tau */
tau_3(ai,aj,ak,lambda=Pi/2,verbose=0) = {
  my(zk,xItau);
  /* inverse ai->-1, aj->1 */
  zk = (2*ak-ai-aj)/(aj-ai);
  if(verbose,print("zk=",zk));
  xItau = asinh(atanh(zk)/lambda);
  if(verbose,print("x+I*tau=",xItau));
  return(abs(imag(xItau)));
}
tau_edge(A,i,j,lambda=Pi/2) = {
 tauij = 4; /* > Pi/2 */
  for(k=1,#A,
    if(k==i||k==j,next());
      tauij = min(tauij,tau_3(A[i],A[j],A[k]));
      );
   return(tauij);
}
tau_tree(A,edges) = {
  my(tau = 4);
  for(k=1,#edges,
    tau = min(tau,tau_edge(A,edges[k][1],edges[k][2]));
  );
  return(tau);
}

/* bound M1 on g */
dist_1(p) = {
  my(xp = abs(real(p)));
  if(xp>1,abs(xp-1+I*imag(p)),abs(imag(p)));
}
bound_M1(Aprim) = {
  1/sqrt(prod(k=1,#Aprim,dist_1(Aprim[k])));
}
t_opt(D,M1,lambda=Pi/2)=asinh((D+log(4*M1))/lambda);

/* bound M2 on Δ_τ */
dist_thsh(p,tau,lambda=Pi/2) = {
  \\ take everything in the first quatran
  p = abs(real(p))+I*abs(imag(p));
  x0 = 0; x1 = acosh(Pi/(2*lambda*sin(tau)));\\ st. λch(x)sin(τ)=Pi/2
  \\ solve orthogonality (p-φ(x)).φ'(x) = 0,
  \\ but with arguments
  my (phi(x)=tanh(lambda*sinh(x+I*tau)));
  (argphip(x) = arg(cosh(x+I*tau))-2*arg(cosh(lambda*sinh(x+I*tau))));
  (argphi(x)=arg(p-phi(x)));
  x = solve(x=x0,x1, argphi(x)-argphip(x)-Pi/2);
  abs(p-phi(x));
}
bound_M2(Aprim,tau) = {
  1/sqrt(prod(k=1,#Aprim,dist_thsh(Aprim[k],tau)));
}
/* int_R lambda*cht/ch(lambda*sh t) < 4, 1/sqrt(g(phi)) < M */
/* majoration de \int λ\ch(t)/\ch(λ\sh(t+iτ))^(2α) */
phi_bound(tau,lambda=Pi/2,alpha=1/2) = {
  \\ fixed value of Ytau
  Ytau = sqrt(Pi/2*lambda*sin(tau));
  Xtau = sqrt(Ytau^2-lambda^2*sin(tau)^2)/tan(tau);
  return((1+2*alpha)/(2*alpha*cos(Ytau)^(2*alpha))+1/(2*alpha*sinh(Xtau)^(2*alpha)));
  /* it is not very bad compared to the better value below */
  \\Ytau=solve(y=lambda*sin(tau),Pi/2.001,2/cos(y)-log(exp(sqrt(y^2-lambda^2*sin(tau)^2)/tan(tau))-1));
  \\return(4/cos(Ytau));
}
h_opt(D,tau,M2,lambda=Pi/2)=2*Pi*tau/log(1+2*M2*phi_bound(tau,lambda,0.5)*exp(D));

/* make Aprim vectors : remove one or two roots */
make_Aprim(A,ia,ib=0,p) = {
  if(ib, \\ between a(-> -1) and b (-> 1) Weierstrass points
    my(Aprim=vector(#A-2));
    my(a = A[ia]); my(b=A[ib]); my(j=0);
    for(i=1,#A,i==ia||i==ib||Aprim[j++]=(2*A[i]-a-b)/(b-a));
    return(Aprim);
   , \\ a Weierstrass point -> -1 ; p -> 1
    my(Aprim=vector(#A-1));
    my(a=A[ia]);my(j=0);
    for(i=1,#A,i==ia||Aprim[j++]=(2*A[i]-a-p)/(p-a));
    return(Aprim);
    );
}

integration_parameters(A,tree,tau,D,provenbounds=0) = {
  if(provenbounds,
    \\ prove all computations
    tau *= .95; \\ do it precisely
    M1 = M2 = 0;
    for(i=1,#tree,
      i1 = tree[i][1]; i2 = tree[i][2];  
      my(Aprim = make_Aprim(A,i1,i2));
      M1 = max(M1,bound_M1(Aprim));
      M2 = max(M2,bound_M2(Aprim,tau));
      );
    msgtimer("bounds");
    ,
    \\ otherwise this should work, not proven
    M1 = M2 = 1; \\ not true
    );
  h = h_opt(D,tau,M2);
  npoints = ceil(t_opt(D,M1)/h);
  return([h,npoints]);
}

/**
precompute integration points 
   x = [ tanh(λ sinh(kh)) ]
  dx = [ λ cosh(kh)/cosh(λsinh(kh))
suitable for integration on [-1,1]
  _with simplification by a factor cosh(λsinh)_.
*/
integration_points_thsh(h,npoints,lambda=Pi/2) = {
	my(k,eh,eh_inv,res,ekh,ekh_inv,sh,chi,thsh,supres);
	eh = exp(h);
	eh_inv = 1/eh;
	ekh = ekh_inv = 1;
	res = vector(npoints); 
	for(k=1,npoints,
			ekh *= eh;
			ekh_inv *= eh_inv;
      sh = (ekh-ekh_inv)/2;
      ch2 = (ekh+ekh_inv); /* 2*cosh(kh) */
      esh = exp(lambda*sh);
      esh_inv = 1/esh;
      chsh2_i = 1/(esh+esh_inv); \\ 1/(2*cosh(lambda*sinh(kh)))
      shsh2 = esh-esh_inv; \\ 2*sinh(lambda*sinh(kh))
      thsh = shsh2*chsh2_i;
      res[k] = [thsh,ch2*chsh2_i];
     );
  return([h*lambda,res]);
}

/**
  -- The problem of the square root --
 */

/* transform first the problem to [-1,1] where the determination issue is
 * easy */

/* square root with standard determination for each monomial part */
/* use standard transform to express cos(θ/2) and sin(θ/2), see
   "How to Find the Square Root of a Complex Number" by Stanley Rabinowitz
*/

/* sign of the imaginary part with positive zero */
isign(x) = { my(s);s=sign(imag(x));if(s,s,1);}

/* we want to compute sqrt(z-a1)*...sqrt(z-ai)
   with only one square root sqrt((z-a1)*...(z-ai)),
   so we need to track the number of loops around 0
   while multiplying the (z-ai).
   The product flag is for testing purposes.
   TODO : choose square root method
 */
sqrt_affinereduction(Aprim,z,product=0) = {
  my(z1,i1,z2,i2,sgn);
  product=1; \\ FIXME should choose it only once
  \\ if we want the long method, or there is only one root
    if(product||#Aprim==1,
        return(prod(i=1,#Aprim,
            if(real(Aprim[i])>0,I*sqrt(Aprim[i]-z),sqrt(z-Aprim[i]))
            ))
      );
  \\ otherwise
    sgn = 1; \\ track the number of loops
    \\if(real(Aprim[1])>0,z1=Aprim[1]-z;sgn*=-I,z1=z-Aprim[1]);
  z1 = z-Aprim[1];
  /* if the point is on the right, take
     non-standard cut */
  if(imag(z1)==0&&real(z1)<0,z1=-z1;sgn*=-I);
  i1 = isign(z1);
  for(i=2,#Aprim,
      \\if(real(Aprim[i])>0,z2=Aprim[i]-z;sgn*=-I,z2=z-Aprim[i]);
      z2 = z-Aprim[i];
      if(imag(z2)==0&&real(z2)<0,z2=-z2;sgn *=-I);
      i2 = isign(z2);
      z1 *= z2;
      if(i1==i2,
        i1 = isign(z1);
        if(i1!=i2,sgn*=-1),
        i1 = isign(z1)
        );
     );
  return(sqrt(z1)*sgn);
}


/* compute the column of period integrals over [A[i1],A[i2]]
   [ int x^k/y dx ]
*/

/* integral over a loop enclosing A[i1] and A[i2], up to a sign */
/* assume #A>=3 */
int_periods_affinereduction(C,i1,i2,decomp=0) = {
  my(A,a,b,dec,int_points,geom_factor,decprim,Aprim,res,x,dx,tp,tm,j);
  A = C[iRoots];
  int_points = C[iIntPoints];
  a = A[i1]; b = A[i2];
  \\ always integrate from left to right : this ensures \sqrt{b-a} to be
  \\ well defined
  if(real(a)>real(b),print("bad order of edge"));
  dec = (a+b)/2;
  geom_factor = (b-a)/2;
  decprim = (b+a)/(b-a);
  /* affine transformation to [-1,1] */
  Aprim = make_Aprim(A,i1,i2);
  msgdebug([a,b], "=========== de ==========\n    a,b =", 2);
  /* numerical integration with precomputed DE weights */
  \\ init on x = 0, dx=1
  my(g=C[iGenus]);
  res = vector(g);
  res[1] = 1/sqrt_affinereduction(Aprim,0,decomp);
  for(j=1,g-1,res[j+1] = res[j] * decprim);
  \\ integrate on precomputed points
  for(i=1,#int_points,
      x = int_points[i][1];
      dx = int_points[i][2];
      tp  = dx/sqrt_affinereduction(Aprim,+x,decomp);
      tm = dx/sqrt_affinereduction(Aprim,-x,decomp);
        res[1] += (tp+tm);
      for(j=2,g,
        tp *= (decprim+x);
        tm *= (decprim-x);
        res[j] += (tp+tm);
      );
     );
   res *= C[iIntFactor];
   \\ factor due to the change of variable
   msgdebug(res, "    shifted =", 2);
   fact = I/sqrt(geom_factor)^#Aprim;
   res[1] *= fact;
   for(j=2,g,
     fact *= geom_factor;
     res[j] *= fact;
     );
   msgdebug(res, "    factors =", 3);
   \\ to track signs, return also the values of the sqrt chosen at end points
   vplus = sqrt_affinereduction(Aprim,1,decomp)*sqrt(geom_factor)^#Aprim/geom_factor;
   vmoins = sqrt_affinereduction(Aprim,-1,decomp)*sqrt(geom_factor)^#Aprim/geom_factor;
   return([res,vmoins,vplus]);
}

/* integrate between a and p */
newint_AbelJacobi_affinereduction(C,i1,p,decomp=0) = {
  my(A=C[iRoots]);
  my(g=C[iGenus]);
  return(intnum(x=[A[i1],0.5],p,vector(g,k,z=prod(i=1,#A,sqrt(x-A[i]));if(z,x^(k-1)/z,1))));
}
/** partial integration from P to Pi.
Assume that no branch point is on [P Pi]
*/
int_AbelJacobi_affinereduction(C,i1,p,y_p,decomp=0) = {
  my(A=C[iRoots]);
  my(int_points=C[iIntPoints]);
  my(int_factor=C[iIntFactor]);
  a=A[i1];
  dec = (a+p)/2;
  geom_factor = (p-a)/2;
  decprim = (p+a)/(p-a);
  /* affine transformation to [-1,1] */
  Aprim = make_Aprim(A,i1,,p);
  /* numerical integration with precomputed DE weights */
  \\ init on x = 0, dx=1
  my(g=C[iGenus]);
  res = vector(g);
  res[1] = 1/sqrt_affinereduction(Aprim,0,decomp);
  for(j=1,g-1,res[j+1] = res[j] * decprim);
  \\ integrate on precomputed points
    for(i=1,#int_points,
        x = int_points[i][1];
        dx = int_points[i][2];
        tp  = dx/sqrt_affinereduction(Aprim,+x,decomp)*sqrt(1-x);
        tm = dx/sqrt_affinereduction(Aprim,-x,decomp)*sqrt(x+1);
        res[1] += (tp+tm);
        for(j=2,g,
          tp *= (decprim+x);
          tm *= (decprim-x);
          res[j] += (tp+tm);
          );
       );
  \\ factor due to the change of variable
    fact = -int_factor/sqrt(geom_factor)^(#Aprim-1);
  res[1] *= fact;
  for(j=2,g,
      fact *= geom_factor;
      res[j] *= fact;
     );
  /* right sheet ? */
  y_p1 = sqrt((p-a)/2)^(#A)*sqrt_affinereduction(Aprim,1,decomp);
  sgn = sign(real(y_p/y_p1));
  return(sgn*res);
  \\x_p = prod(i=1,#A,p-A[i]);
  \\ FIXME : remove below
  
  /* find the right sheet : the result has to be the same as if computed using
 the point P0 */
  /* FIXME : does not work if P0P passes through some Pi */
  \\ compute using [Pi p]
  fiP = sqrt((p-a))^#A*sqrt_affinereduction(Aprim,1,decomp);
  \\ compute using [P0 p] ...
  P0 = A[C[iP0]];
  f0P = sqrt((p-P0))^#A*sqrt_affinereduction(make_Aprim(A,C[iP0],,p),1,decomp);
  sgn = sign(real(fiP/f0P));
  msgdebug(Strprintf("sheets : %+1.3f",real(fiP/f0P)));
  return(sgn*res);
}
/* find point on a conventional upper sheet given its x-coordinate */
hcUpperSheet(C,x,flag=0) = {
  my(A = C[iRoots]);
  if(!flag,
    \\ standard cut at each branch point
    return([x,prod(i=1,#A,sqrt(x-A[i]))])
  ,if(flag==1,
    \\ global standard cut
    return([x,sqrt(prod(i=1,#A,x-A[i]))])
   ));
  }
/* return the intersection number of oriented [ab] and [cd] */
intersection_intervals(A,ia,ib,ic,id) = {
  \\ check no double point
  \\if((a-b)*(a-c)*(a-d)*(b-c)*(b-d)*(c-d)==0,return(0));
  a = A[ia]; b = A[ib]; c = A[ic]; d = A[id];
  cprim = (2*c-a-b)/(b-a);
  dprim = (2*d-a-b)/(b-a);
  imc = imag(cprim);
  imd = imag(dprim);
  if(sign(imc)*sign(imd)==1,return(0)); \\ on the same side
  /* p the intersection */
  xp = imag(conj(cprim)*(dprim-cprim));
  fp = imd-imc;
  if(abs(xp)>abs(fp),return(0)); \\ intersection at xp/fp in [-1,1] ?
  \\ result
  return(sign(fp));
}

/* return the intersection number between the periods
   I_α and I_β computed with standard formula */
/* this works except with aligned points ... */
intersection(A,ia,ib,ic,id) = {
  \\ do we have a common point ?
  if(ia==ib||ic==id,return(0)); \\ bad entry
  if((ia==ic&&ib==id)||(ia==id&&ib==ic),return(0)); \\ self intersection
  if(ia==ic,return(+intersection_abad(A,ia,ib,ic,id)));
  if(ib==ic,return(+intersection_abbd(A,ia,ib,ic,id)));
  if(ia==id,return(-intersection_abbd(A,ic,id,ia,ib)));
  if(ib==id,return(+intersection_abcb(A,ia,ib,ic,id)));
  \\ transverse case, cannot occur for a maximum tree
  return(intersection_inner(A,ia,ib,ic,id));
}

/* end intersection I[ab].I[bd] */
intersection_abbd(A,ia,ib,ic,id) = {
  k = #A;
  a = A[ia]; b = A[ib]; d = A[id];
  tau = arg((d-b)/(b-a)); 
  /* \frac{\sqrt{d-b}^{k}\prod_{λ\neq b,d}\sqrt{x_β(b)-ψ_β(λ)}} */
  \\fbd = sqrt((d-b))^k
  \\*prod(i=1,k,if(i==ib||i==id,1,sqrt(-1-(2*A[i]-b-d)/(d-b))));
  Aprim = make_Aprim(A,ib,id);
  fbd = sqrt((d-b))^k*sqrt_affinereduction(Aprim,-1);
  \\fab = sqrt((b-a))^k
  \\*prod(i=1,k,if(i==ib||i==ia,1,sqrt(+1-(2*A[i]-b-a)/(b-a))));
  Aprim = make_Aprim(A,ia,ib);
  fab = sqrt((b-a))^k*sqrt_affinereduction(Aprim,1);
  if(default(debug)>=1,
      endarg = arg(fab/fbd);
      arg1 = -(Pi+tau)/2;
      arg2 = (Pi-tau)/2;
        printf("intersection_abbd %i->%i->%i\n",ia,ib,id);
        print("  fab=",fab);
        print("  fbd=",fbd);
        print("  tau=",tau);
        printf("  arg[fab/fbd]=%1.4f should be %1.4f or %1.4f\n",endarg,arg1, arg2);
        print("   roots ab ", Aprim);
      if(abs(endarg-arg1)>0.0001&&abs(endarg-arg2)>0.0001,
        printf("tau/2=%1.3f ; 2*tau=%1.3f ; (Pi-tau)/2=%1.3f ; Pi-2*tau=%1.3f ; \n",
          tau/2,2*tau,(Pi-tau)/2,Pi-2*tau);
        );
    );
  if(imag(fbd/fab)<0,return(-1),return(1));
}
intersection_abcb(A,ia,ib,ic,id) = {
  k = #A;
  a = A[ia]; b = A[ib]; c = A[ic];
  tau = arg((b-c)/(b-a)); 
  \\ \frac{\sqrt{d-b}^{k}\prod_{λ\neq b,d}\sqrt{x_β(b)-ψ_β(λ)}}
  \\fcb = sqrt((b-c))^k
    \\*prod(i=1,k,if(i==ic||i==ib,1,sqrt(+1-(2*A[i]-c-b)/(b-c))));
  Aprim = make_Aprim(A,ic,ib);
  fcb = sqrt((b-c))^k * sqrt_affinereduction(Aprim,+1);
  \\fab = sqrt((b-a))^k
    \\*prod(i=1,k,if(i==ib||i==ia,1,sqrt(+1-(2*A[i]-b-a)/(b-a))));
  Aprim = make_Aprim(A,ia,ib);
  fab = sqrt((b-a))^k*sqrt_affinereduction(Aprim,+1);
  if(default(debug)>=1,
      endarg = arg(fab/fcb);
      \\arg1 = tau/2;
      arg1 = if(tau>0,Pi-tau/2,-Pi-tau/2);
      \\arg2 = Pi-tau/2;
      arg2 = -tau/2;
        print("intersection_abcb");
        print("  tau=",tau);
        printf("  arg[fab/fcb]=%1.4f should be %1.4f or %1.4f\n",endarg,arg1, arg2);
      if(abs(endarg-arg1)>0.0001&&abs(endarg-arg2)>0.0001,
        printf("tau/2=%1.3f ; 2*tau=%1.3f ; (Pi-tau)/2=%1.3f ; Pi-2*tau=%1.3f ; Pi-tau/2=%1.3f ; Pi+tau/2=%1.3f \n", tau/2,2*tau,(Pi-tau)/2,Pi-2*tau,Pi-tau/2,Pi+tau/2);
        );
    );
  if(imag(fab/fcb)<0,return(-1),return(1));
}
intersection_abad(A,ia,ib,ic,id) = {
  k = #A;
  a = A[ia]; b = A[ib]; d = A[id];
  tau = arg((d-a)/(b-a)); 
  \\ \frac{\sqrt{d-b}^{k}\prod_{λ\neq b,d}\sqrt{x_β(b)-ψ_β(λ)}}
  \\fad = sqrt((d-a))^k
    \\*prod(i=1,k,if(i==ia||i==id,1,sqrt(-1-(2*A[i]-a-d)/(d-a))));
  Aprim = make_Aprim(A,ia,id);
  fad = sqrt((d-a))^k*sqrt_affinereduction(Aprim,-1);
  \\fab = sqrt((b-a))^k
    \\*prod(i=1,k,if(i==ib||i==ia,1,sqrt(-1-(2*A[i]-b-a)/(b-a))));
  Aprim = make_Aprim(A,ia,ib);
  j=0;for(i=1,k,if(i==ia||i==ib,,Aprim[j++]=(2*A[i]-a-b)/(b-a)));
  fab = sqrt((b-a))^k*sqrt_affinereduction(Aprim,-1);
  if(default(debug)>=1,
      endarg = arg(fab/fad);
      arg1 = -tau/2;
      arg2 = if(tau>0,Pi-tau/2,-Pi-tau/2);
        print("intersection_abad");
        print("  tau=",tau);
        printf("  arg[fab/fad]=%1.4f should be %1.4f or %1.4f\n",endarg,arg1, arg2);
      if(abs(endarg-arg1)>0.0001&&abs(endarg-arg2)>0.0001,
        printf("tau/2=%1.3f ; 2*tau=%1.3f ; (Pi-tau)/2=%1.3f ; Pi-2*tau=%1.3f ; Pi-tau/2=%1.3f ; Pi+tau/2=%1.3f \n", tau/2,2*tau,(Pi-tau)/2,Pi-2*tau,Pi-tau/2,Pi+tau/2);
        );
    );
  if(imag(fab/fad)<0,return(-1),return(1));
}


/* inner intersection I[ab].I[cd] */
/* assume different end points */
intersection_inner(A,ia,ib,ic,id) = {
  k = #A;
  a = A[ia]; b = A[ib]; c = A[ic]; d = A[id];
  /* first rule out null intersection */
  cprim = (2*c-a-b)/(b-a);
  dprim = (2*d-a-b)/(b-a);
  imc = imag(cprim);
  imd = imag(dprim);
  if(sign(imc)*sign(imd)==1,return(0)); \\ on the same side
    /* p the intersection */
    xp = imag(conj(cprim)*(dprim-cprim));
  fp = imd-imc;
  if(abs(xp)>=abs(fp),return(0)); \\ discard if xp not in ]-1,1[
    xpab = xp/fp;
  \\ compare the sheets at the intersection p */
    p = (b-a)/2*xpab+(b+a)/2;
  xpcd = (2*p-c-d)/(d-c); \\ should be in ]-1,1[
    /* value 1/y(p) on [ab] */
    fpab = sqrt(b-a)^#A * prod(i=1,#A,if(i==ia||i==ib,1,sqrt(xpab-(2*A[i]-a-b)/(b-a))));
  /* value 1/y(p) on [cd] */
  fpcd = sqrt(d-c)^#A * prod(i=1,#A,if(i==ic||i==id,1,sqrt(xpcd-(2*A[i]-c-d)/(d-c))));
  if(default(debug)>2,
      printf("intersection_inner "
        "[%1.2f%+1.2f*I,%1.2f%+1.2f*I].[%1.2f%+1.2f*I,%1.2f%+1.2f*I]\n",
      real(a),imag(a),real(b),imag(b),real(c),imag(c),real(d),imag(d));
      \\ check everything is correct
      /* xpcd */
      aprim = (2*a-c-d)/(d-c);
      bprim = (2*b-c-d)/(d-c);
      xpcd2 = imag(conj(aprim)*(bprim-aprim))/imag(bprim-aprim);
      print("  xpcd=",xpcd,"=",xpcd2);
      /* quotient and cosh quotient */
      tpab = atanh(xpab);
      tpcd = atanh(xpcd);
      print("  fpab/fpcd=",fpab/fpcd);
      print("  cosh[tpab]/cosh[tpcd]=",cosh(tpab)/cosh(tpcd));
      refab = I*cosh(tpab)/fpab*2^(#A/2); refp = 1/prod(i=1,#A,sqrt(p-A[i]));
      printf("  and fpab=%1.3f%+1.3fI should be +/-f[p]=%1.3f%+1.3fI\n",
        real(refab),imag(refab),real(refp),imag(refp));
    );
  \\ result : sign(fp) gives intersection number [ab].[cd]
    \\ if sgn(fpab/fpcd)=+1
    return(2*sign(fp)*sign(real(fpab/fpcd)));
}

/* *************************************************************************
 *  period matrix with maximal-flow spanning tree
 * ************************************************************************* */

/* step 1 : set of maximum spanning-tree edges and its flow */
\\max_spanning(A,nedges=0)= {
\\  my(n,tau_v,z,taken,tree,taumin);
\\  n = #A;
\\  if(!nedges,nedges=n-1); /* spanning tree */
\\  /* start from A[1] */
\\  con = vector(n,i,0); /* connected vectices */
\\  con[1] = 1; nc = 1;
\\  for(nc = 2, n-1,
\\     
\\  );
\\  /* if nedges just remove one */
\\}

 
\\ FIXME : very bad algorithm for large genus
max_spanning(A,nedges=0)= {
  my(n,tau_v,z,taken,tree,taumin);
  n = #A;
  if(!nedges,nedges=n-1); /* spanning tree */
  tau_v=vector(n*(n-1)/2);
  /* build list of edges ordered by tau */
  z=0;
  for(i=1,n,
    for(j=i+1,n,
      tau_v[z++] = [tau_edge(A,i,j),i,j];
      \\ order edge so that left to right
      if(real(A[i])>real(A[j]),tau_v[z][2]=j;tau_v[z][3]=i);
      );
    );
  tau_v = vecsort(tau_v,1); \\ sort according to tau
  taken = vector(n,i,0); \\ mark vertex if in tree
  tree = vector(nedges);
  /* this is not optimal... too many loops */
  z=n*(n-1)/2; \\select edge with maximum tau
  taken[tau_v[z][2]]=1; \\ z taken at first loop
  taumin = tau_v[z][1];
  for(k=1,nedges,
    z=n*(n-1)/2;
    /* consider next edge with exactly one vertex taken (no cycle) */
    while(taken[tau_v[z][2]]+taken[tau_v[z][3]]!=1,z--);
    /* this is the best edge allowed */
    i = tau_v[z][2]; j = tau_v[z][3];
    tree[k] = [i,j];
    taumin = min(taumin,tau_v[z][1]);
    taken[i] = taken[j] = 1;
    );
  return([tree,taumin]);
}

/* step 2 : intersection matrix of selected edges */
intersection_spanning(A,tree) = {
  n = #tree;
  res = matrix(n,n);
  for(i=1,#tree,
    for(j=i+1,#tree,
    \\print("intersection ",i,".",j,":");
    s = intersection(A,tree[i][1],tree[i][2],tree[j][1],tree[j][2]);
    res[i,j] =  s;
    res[j,i] = -s;
    );
  );
  return(res);
}

/* step 3 : period matrix of edges + intersection matrix */
/* suppose integration points precomputed */
periods_spanning(C) = {
  my(g = C[iGenus]);
  my(tree = C[iTree]);
  my(A = C[iRoots]);
  /* nb : should give tau since already computed.. */
  if(keependpts=1,
  Mat(vector(#tree,k,int_periods_affinereduction(C,tree[k][1],tree[k][2])[1]~));
  ,
  res = vector(#tree);
  endpointsvalues=vector(#tree);
  for(k=1,#tree,
    tmp = int_periods_affinereduction(A,tree[k][1],tree[k][2]);
    /* should use final values to find intersection if sure that no
     * crossing edges TODO */
    res[k] = tmp[1]~;
    endpointsvalues[k] = [tmp[2],tmp[3]];
  );
  Mat(res);
  );
}

/* step 4 : retrieve symplectic basis */

/* f[i] <- u*f[i]+v*f[j] */
\\column_operation(i,j,u,v) = {
\\  }

/* assume M is invertible, antisymetric and of even dimension 2g */
symplectic_reduction(M) = {
  g=matsize(M)[1];
  if(g%2,return(0),g/=2);
  N = Vec(M); \\ columns matrices
  P = Vec(matid(2*g));
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
  /* and finally canoniq basis */
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

/* choose point such that tau is maximal */
choose_integration_path(C,p) = {
  my(besttau = 0); /* > Pi/2 */
  my(tau = C[iIntTau]);
  my(A = C[iRoots]);
  for(i=1,#A, \\ choose A[i]
      tauk = 4; /* > Pi/2 */
      for(k=1,#A,
        if(k==i,next());
        tauk = min(tauk,tau_3(p,A[i],A[k]));
        );
      if(tauk>besttau, besti=i; besttau = tauk);
      if(besttau>tau,
        msgdebug(Str("choose point a[",besti,"]=",A[besti]));
        return([besti,besttau])
        ); \\ OK with precomputed points
      );
   \\ ah, need to subdivide
   \\ TODO
   msgdebug("bad accuracy : needs subdivision");
   return([besti,besttau]);
}
hcModABlattice(C,z) = {
  my(g=C[iGenus]);
  my(Rproj=C[iReduce]);
  Rz = concat(real(z),imag(z));
  centerz=frac(Rproj*Rz);
  return(centerz);
}

/** compute Abel-Jacobi map of the point P = (x,y)


*/
/* give image in basis (A,B) with standard periods */
/* give image in [Id+tau] \R^{2g}/\Z^{2g} */
hcAbelJacobi(C,p) = {
  /* is p one of the roots ? */
  my(A=C[iRoots]);
  my(res);
  if(type(p)=="t_VEC",
    x_p = p[1];
    y_p = p[2]
  ,
    x_p = p;
    y_p = hcUpperSheet(C,p)[2];
    );
  for(i=1,#A,if(x_p==A[i],
        /* it is equal to half the period P0-P */
        return(hcWeierstrassInt(C)[,i]);
        );
     );
  /* otherwise perform integration */
  my(tmp = choose_integration_path(C,x_p));
  res = int_AbelJacobi_affinereduction(C,tmp[1],x_p,y_p);
  /* check the sheet */
  res = hcModABlattice(C,res~);
  return(frac(res+hcWeierstrassInt(C)[,tmp[1]]));
  \\  /* choose closest point */
  \\  dmin=abs(p-A[1]);
  \\  imin=1;
  \\  for(i=2,#A,if(u=abs(p-A[i])<dmin,dmin=u;imin=i;break()));
  \\  a=A[imin];
  \\  taumin=3;
  \\  for(i=1,#A,if(i!=imin,taumin=min(taumin,tau_3(p,a,A[i]))));
  \\  /* the path should be optimized and validated */
  \\  /* integrate */
  \\  int_AbelJacobi_affinereduction(C,imin,p);
}
/* find the integrals between P0 and Pi using P */ 
find_hcWeierstrass(C,p) = {
  my(g=C[iGenus]);
  my(P0=C[iP0]);
  P0P = hcModABlattice(C,int_AbelJacobi_affinereduction(C,P0,p)~);
  round(2*Mat(vector(2*g+1,k,
  hcModABlattice(C,int_AbelJacobi_affinereduction(C,k,p)~)
    -P0P)))/2;
} 

/* *************************************************************************
 *  put everything together
 * ************************************************************************* */

/** Remark on matrices conventions : {{{
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
}}}
*/

/** image of Weierstrass points under Abel-Jacobi map
*/
AJWeierstrass(C,Cbasis=0) = {
  /* we take one Weierstrass point as origin and derive the image of all other
following the path given by the spanning tree.
 */
  my(A=C[iRoots]);
  my(g=C[iGenus]);
  /* matrix of paths P0Pi in basis C_1,.. C_2g */
  my(P0Pi = matrix(2*g,#A,i,j,0));
  my(taken=vector(#A,k,0));
  my(tree=C[iTree]);
  /* first take the first point as P0 */
  my(P0 = tree[1][1]); taken[P0] =1;
  for(i=1,#tree,
      my(Pstart,Pend);
      /* fill the new point */
      if(taken[tree[i][1]],
      Pstart = tree[i][1];
      Pend   = tree[i][2];
      taken[tree[i][2]]=1
      ,
      Pstart = tree[i][2];
      Pend   = tree[i][1];
      taken[tree[i][1]]=1
      );
      P0Pi[,Pend] = P0Pi[,Pstart];
      P0Pi[i,Pend] = 1/2
     );
  /* substract the real P0 column from all */
  P0 = C[iP0];
  \\P0P0 = P0Pi[,P0]; /* also works */
  P0P0 = vector(2*g,k,P0Pi[k,P0])~;
  for(j=1,#A,
      P0Pi[,j] = P0Pi[,j] - P0P0;
     );
  if(Cbasis,return(Mat(P0Pi)));
  /* and in basis A,B */
  my(ABtoC=C[iABtoC]);
  return(frac(ABtoC^(-1)*Mat(P0Pi)));
}

/** compute Abel-Jacobi map of Weierstrass points
   in the standard Mumford AB basis */
AJWeierstrass_Mumford(C) = {
  my(A=C[iRoots]);
  my(g=C[iGenus]);
  my(P0=C[iP0]);
  /*  second method : same assumption as above; by definition
      we know the expression of loops A_i and B_i in terms of P_i
      and just invert the matrix.
   */
  \\for(i=1,g,
  \\  ABofPi[g+i,2*i-1]=2; ABofPi[g+i,2*i]=2;
  \\  for(j=2*i,2*g+1,ABofPi[i,j]+=2)
  /* old without inversion */
  \\ A,B is the standard symplectic basis so we know
  \\ the image of each difference of points
  P0Pi = matrix(2*g,#A);
  /* FIXME : A_i and B_i have been exchanged */
  for(i=1,g,
    \\ 1/2 A_i = a_{2i-1}+a_{2i} (mod A_i)
    P0Pi[g+i,2*i-1]=1/2; P0Pi[g+i,2*i]=1/2;
    \\ 1/2 B_i = a_{2*i}+a_{2*i+1}+..a{2*g+1}
    for(j=2*i,2*g+1,P0Pi[i,j]+=1/2)
  );
  \\ and substract P0
  P0P0 = P0Pi[,P0];
  for(j=1,#A,
      P0Pi[,j] = P0Pi[,j] - P0P0;
     );
  return(frac(P0Pi));
}
/** homology base change from the computed symplectic basis into Mumford A,B
 basis defined with the ordered set of roots. 
 The change of basis is determined using the image of Weierstrass points

   P1   P2g+1        C1     C2g        P1   P2g+1
  +----------+      +----------+      +----------+
A1|          |    A1|          |   C1 |          |
  |          | =    |          | *    |          |
Bg|          |    Bg|          |   C2g|          |

     W_in_AB =      M_CtoAB      *    W_in_C   mod \Z^2g
ie we have M_CtoAB in \F2 and we have to lift it to \Z
such that it is in Sp(g).
 */
ABtoMumfordAB_mod(C) = {
  /* assume P1 is the origin */
  P0 = C[iP0];
  if(P0==1,
    strcol = Str("2..")
    ,if(P0==2*g+1,
    strcol = Str("1..-2")
    ,
    strcol = Str("1..",P0-1,",",P0+1,"..",2*g+1)
    ));
    strlin = "1..";
  WinC_AB = vecextract(2*C[iAJRoots],strlin,strcol);
  WinMumfordAB = vecextract(2*AJWeierstrass_Mumford(C),strlin,strcol);
  return((WinMumfordAB*WinC_AB^(-1))%2)
}
/* the real matrix */
ABtoMumfordAB(C) = {
  return(MumfordABtoC(C)^(-1)*C[iABtoC]);
}
/** express C basis in loops A_i, B'_i, modulo 2\Z... */
CtoABp(C) = {
  g = C[iGenus];
  tree = C[iTree];
  M = matrix(2*g,2*g);
  for(j=1,2*g,
    if(tree[j][1]<tree[j][2],s=tree[j][1];e=tree[j][2],
                             e=tree[j][1];s=tree[j][2]);
    for(iA = floor((s+2)/2),floor(e/2),
       M[iA,j] = 1
       );
    for(iBp = floor((s+1)/2), floor((e-1)/2),
       M[g+iBp,j] = 1
       );
    );
  return(M);
}
CtoAB(C) = {
  return(ABptoAB(C[iGenus]) * CtoABp(C));
}
 
/** intersection between A_i, B'_i and C_j as I would compute them */
int_ABp_C(C) = {
  g = C[iGenus];
  tree = C[iTree];
  A = C[iRoots];
  M = matrix(2*g,2*g);
  for(j=1,2*g,
    s = tree[j][1];
    e = tree[j][2];
    for(i=1,g,
      M[i,j]   = intersection(A,2*i-1,2*i,s,e);
      M[g+i,j] = intersection(A,2*i,2*i+1,s,e);
       );
     );
  return(M);
}
/** intersection matrix of loops A_i, B_i' as I _would_ compute them */
int_ABp(C) = {
  my(g = C[iGenus]);
  my(A = C[iRoots]);
  my(M = matrix(2*g,2*g));
  for(i=1,g,
    for(j=1,g,
      /* A_i.A_j */
      M[  i,  j] = intersection(A,2*i-1,2*i,2*j-1,2*j);
      /* A_i.B_j */
      M[  i,g+j] = intersection(A,2*i-1,2*i,2*j,2*j+1);
      /* B_i.A_j */
      M[g+i,  j] = intersection(A,2*i,2*i+1,2*j-1,2*j);
      /* B_i.B_j */
      M[g+i,g+j] = intersection(A,2*i,2*i+1,2*j,2*j+1);
   );
  );
return(M);
}
/* if all signs were compatible the above would be equal to the following. It
 * must be equal mod 2\Z */
int_ABp_ref(x) = {
  my(g); if(type(x)=="t_VEC", g = C[iGenus], g=x);
  my(P = ABptoAB(g));
  return(P*matJ(g)*mattranspose(P));
  return(mattranspose(P)*matJ(g)*P);
}

ABtoABp(g) = {
  my(P = matrix(2*g,2*g));
  for(j=1,g,
    /* A_j = A_j */
    P[j,j] = 1;
    /* B_j = B_j'+...B_g' */
    for(i=j,g,P[g+i,g+j]=1)
    );
  return(P);
}
ABptoAB(g) = {
    /* A_i = A_i */
    /* B_i' = B_i-B_{i+1} */
  my(P = matid(2*g));
  for(i=1,g-1,
     P[g+i,g+i+1] = -1
     );
  return(P);
}
/** test all signs on C until the A,B basis is symplectic */
test_all_signs_C(C) = {
  my(g = C[iGenus]);
  my(Ic = C[iIntersection]);
  my(P_C_ABp = CtoABp(C)); 
  my(Iref = int_ABp_ref(g));
  /* vector of signs */
  nsign = trace(mattranspose(P_C_ABp)*P_C_ABp);
  for(i=0,2^nsign-1, \\ Big loop !
     /* fill the signs */
     s=0;j=0;
     while(s<nsign,
          my(l=j \ (2*g)+1);my(c=j%(2*g)+1);
       if(P_C_ABp[l,c]!=0,
          \\print("fill [",l,",",c,"]");
          P_C_ABp[l,c]=(-1)^bittest(i,s);
          s++
          );
        j++
        );
      if(  mattranspose(P_C_ABp) * Iref * P_C_ABp == Ic,
      /* we have found a correct change matrix from C to ABp */
          return(P_C_ABp)
        );
  );
  return("no signs found");
}
/** OK, this time this is correct, however VERY slow algorithm:
The result P is the change of basis between C and Mumford A,B, and
satisfies in particular
 P^t * C[iIntersection] * P = J
Remark : does not work if the Mumford basis intersects itself (ie if it is not
symplectic)
*/
MumfordABtoC(C) = {
  my(P_CtoABp = test_all_signs_C(C));
  return(P_CtoABp^(-1)*ABtoABp(C[iGenus]));
}
 
/* intersection matrix of Mumford loops as I have computed them */
int_AB(C) = {
  my(M = int_ABp(C));
  my(g = C[iGenus]);
  /* base change matrix */
  my(P = ABtoABp(g));
  return(mattranspose(P)*M*P);
  my(P = ABptoAB(g));
  return(P*M*mattranspose(P));
}
isSymplectic(M) = {
 my(g=matsize(M)[1]/2);
 my(a=vecextract(M,2^g-1,2^g-1));
 my(b=vecextract(M,2^g-1,4^g-2^g));
 my(c=vecextract(M,4^g-2^g,2^g-1));
 my(d=vecextract(M,4^g-2^g,4^g-2^g));
 (mattranspose(a)*d-mattranspose(c)*b==matid(g) || print("a^td-c^td!=I"))
 &&
 (mattranspose(a)*c==mattranspose(c)*a || print("a^tc non symmetric"))
 &&
 (mattranspose(b)*d==mattranspose(d)*b || print("b^td non symmetric"))
}

/**
the period matrix is computed in some symplectic homology basis
determined during computations.
*/
hcInit(A,provenbounds=0) = {
  my(hcStruct = vector(iMax));
  my(D,Dcalc,g,tree,tau,h,npoints,IntFactor,IntPoints,AB,Intersection,P);
  my(Omega0,Omega1);
  if(type(A)=="t_POL",A = polroots(A));
  g = floor((#A-1)/2);
  hcStruct[iRoots] = A;
  hcStruct[iGenus] = g;  
  hcStruct[iIntDcalc] = Dcalc = precision(1.);
  /* low precision parameters */
  default(realprecision,19);
  gettime();
  my(tmp = max_spanning(A,2*g));
  msgtimer("tree");
  hcStruct[iTree] =  tree = tmp[1];
  hcStruct[iIntTau] = tau = tmp[2];
  /* prepare integration */
  D = Dcalc*log(10);
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
  hcStruct[iIntegrals] = coh1x_homC;
  /* intersection matrix of basis C_1,.. C_2g */
  IntC = intersection_spanning(A,tree);
  hcStruct[iIntersection] = IntC;
  msgtimer("intersection");
  /* column base change from (A_i,B_i) to C_i */
  /* that is, P^t * Intersection * P = J_g */
  ABtoC = symplectic_reduction(IntC);
  hcStruct[iABtoC] = ABtoC;
  /* change basis */
  /* put everything in the struct */
  coh1x_homAB = coh1x_homC * ABtoC; /* now matrix in (A,B) basis */
  Omega0 = vecextract(coh1x_homAB,Str("1..",g),Str("1..",g));
  Omega1 = vecextract(coh1x_homAB,Str("1..",g),Str(g+1,".."));
  hcStruct[iOmega0] = Omega0; \\ on A basis
  hcStruct[iOmega1] = Omega1; \\ on B basis
  /* with dual basis for H^1 */
  /* WARNING: should one use this van Wamelen convention instead
     of Omega0^-1 Omega1 ?
   */ 
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

/* retrieve values */
hcWeierstrassInt(C) = {
  return(C[iAJRoots])
}
hcBigperiods(C) = {
  return(concat(C[iOmega0],C[iOmega1]));
}
hcSmallperiods(C) = {
  return(C[iTau]);
  return(C[iOmega1]^(-1)*C[iOmega0]);
}
/* intersection matrix of the tree */
hcIntersection(C) = {
  return(C[iIntersection]);
}
/* periods of the tree */
hcTreeperiods(C) = {
  return(hcBigperiods(C)*C[iABtoC]^(-1));
}

/* old conversion functions */
ABtoA_B(AB) = {
  my(A,B);
  g = #AB[,1];
  A = vecextract(AB,Str("1..",g),Str("1..",g));
  B = vecextract(AB,Str("1..",g),Str(g+1,".."));
  return([A,B]);
}
BigtoSmall(AB) = {
  my(tmp = ABtoA_B(AB));
  return(tmp[2]^(-1)*tmp[1]);
}
/**
action of an Moebius transformation on x
*/

homog(H,pts) = {
  res = [];
  for(k=1,#pts,
    my(den = H[2,1]*pts[k]+H[2,2]);
    if(den,res = concat(res,(H[1,1]*pts[k]+H[1,2])/den));
     );
  \\ infinity ?
  if(#pts%2,res=concat(res,0));
  return(res);
  \\return(vector(#pts,k,(H[1,1]*pts[k]+H[1,2])/(H[2,1]*pts[k]+H[2,2])));
}
HS=[0,1;1,0];
HT=[1,1;0,1];
H1=[3,7;2,5];
H2=[8,-3;11,-4];
H3=[13,892;5,343];
H4=[1729,138; 213,17];

/** GL(g) action of a SL(2) moebius transform :
we decompose
(ax+b)/(cx+d) = (b-ad/c)/(cx+d) + a/c
i.e.
[a,b;c,d] = [b-a*d/c,a/c;0,1]*[0,1;1,0]*[c,d;0,1]
*/
/* left matrix corresponding to inversion [0,1;1,0] */
invert_mat(g) = {
  matrix(g,g,i,j,if(i+j==g+1,-1,0));
}
/* left matrix corresponding to homothety ax+b */
homot_mat(g,a,b,degree) = {
  matrix(g,g,i,j,if(i>=j,binomial(i-1,j-1)*a^j*b^(i-j)))/sqrt(a^degree)
}
moebius_matrix_old(m,g,degree) = {
  my(a = m[1,1]);
  my(b = m[1,2]);
  my(c = m[2,1]);
  my(d = m[2,2]);
  if(c,
  return(
homot_mat(g,b-a*d/c,a/c,degree)*invert_mat(g)*homot_mat(g,c,d,degree))
  ,
  if(d,
  return(homot_mat(g,a/d,b/d,degree))
  ,
  return(homot_mat(g,a,b,degree))
  ));
  
}
/* computed all at once */
moebius_matrix(m,g) = {
  my(a = m[1,1]);
  my(b = m[1,2]);
  my(c = m[2,1]);
  my(d = m[2,2]);
/* line k represents (ax+b)^k (cx+d)^(g-k-1) on basis 1,x,...x^{g-1} */
  P1=Pol([a,b]); P2=Pol([c,d]);
  lines = vector(g,kk,
  my(k=kk-1);
  P1^k*P2^(g-1-k)
 );
(a*d-b*c)*
matrix(g,g,k,l,
  polcoeff(lines[k],l-1)
);
}
/** action on the curve */
moebius_action(C,m) = {
  g = C[iGenus];
  A = C[iRoots];
  degree = #A;
  Omega0 = hcBigperiods(C);
  Omega1 = moebius_matrix(m,g)*Omega0;
  Omega2 = hcBigperiods(hcInit(homog(m,A)));
  print("Omega1=",Omega1);
  print("Omega2=",Omega2);
}
