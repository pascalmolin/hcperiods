
/* *************************************************************************
 *  for testing purposes 
 * ************************************************************************* */

/* some hyperelliptic curves */
ptsA = [-3,-1,0,3/2,4]; \\ aligned
ptsB = [-1,2*I,1+I,1-I,-2*I]; \\ hexagon
ptsC = [0,1+I,3*I]; \\ sharp angle

ptsX = [-1,I,-I,1,1+I]; \\ house
ptsX1 = [-1,1,-I,I,1+I]; \\ change order...
ptsX2 = [1+I,-1,1,I,-I];
ptsX3 = [-1,I,1+I,1,-I];
ptsX4 = [-1,1,1+I,-I,I];

ptsEx1 = [0,2-I,-1-I,-3/2+3/2*I,-1+3*I,-3,-3-5/2*I];
ptsEx2 = [5+2*I,5/2+3/2*I,4-I,0,2+5/2*I,-1+3*I,-4+2*I];

ptsHexa = [0,1,1/2+I,-1/2+I,-1,-1/2-I,1/2-I]; \\ hexagon with center
\\ this does a nice tree
{ptsTree= [-101/72 - 175/71*I, 116/51 + 35/19*I, -2/5 + 20/43*I, 2/13 +
1/2*I, 19/23 + 67/79*I, 4/87 - 31/37*I , -11/36 + 82/39*I]}
\\ and same curve with better tau
ptsTreeMoeb=homog([1,1;1,-ptsTree[4]],ptsTree);

\\ different genus
{ptsG7 = [-82/71 - 49/80*I, 409/32 - 7/15*I, 307/108 - 49/47*I, -159/145 -
11/69*I, 13/157 + 175/146*I, 152/51 - 149/189*I, -149/51 - 59/152*I, 142/183 +
5/13*I, -119/97 + 49/24*I, -94/109 + 34/73*I, -70/41 - 9/14*I, 149/169 +
13/120*I, 79/77 + 7/65*I, -5/67 - 149/31*I, 250/167 + 7/23*I]}
{ptsG15 = [-132/37 - 107/40*I, 40/37 - 20/57*I, -562/79 + 93/134*I, 7/152 +
19/21*I, 147/20 - 182/29*I, 206/89 + 167/86*I, -179/147 - 29/150*I, -463/194 -
43/48*I, 655/168 + 31/111*I, -59/160 + 48/115*I, -107/48 + 914/69*I, -178/149 +
304/43*I, 267/109 + 132/103*I, -73/151 + 201/155*I, 154/71 + 20/37*I, 151/196 +
23/75*I, 67/23 - 1433/74*I, 286/163 - 2553/127*I, -134/59 + 1/14*I, 68/97 +
187/64*I, -121/141 + 24/85*I, 1/3 - 1013/71*I, -118/31 + 65/31*I, 13/146 +
79/135*I, -309/110 + 14/51*I, -20/71 + 1/192*I, 529/131 - 12/85*I, 349/144 -
230/81*I, -5/9 - 53/142*I, -275/179 + 141/103*I, 131/80 + 19/36*I]}
{ptsG40 = [25/29 + 20/79*I, 342/157 - 961/162*I, -66/95 - 7/8*I, -348/65 -
458/95*I, -14/29 + 98/33*I, -5/117 - 28/27*I, 34/71 - 134/177*I, 53/163 +
229/16*I, 119/123 - 185/33*I, -5/9 + 311/138*I, 87/59 + 463/89*I, -32/99 -
22/89*I, -160/97 - 161/174*I, 109/185 - 107/23*I, 140/121 + 325/47*I, -21/85 +
311/64*I, -1121/183 - 26/29*I, -160/133 - 57/94*I, 74/11 + 58/63*I, 5/31 -
50/77*I, -60/23 - 257/46*I, 2/9 + 134/133*I, -121/83 - 89/130*I, 21/38 +
55/169*I, -999/125 - 71/30*I, 18/67 + 4/187*I, 153/97 - 91/103*I, 35/96 +
20/71*I, 205/81 + 68/45*I, -100/191 - 23/6*I, 76/181 - 100/127*I, 21/19 -
37/16*I, 292/177 + 133/163*I, -21/8 + 14/41*I, -15/58 + 131/185*I, 621/88 +
460/61*I, -9/184 - 9/89*I, 71/107 - 107/28*I, -39/200 + 23/108*I, 70/33 -
259/127*I, 1/105 + 107/75*I, -233/54 - 9/2*I, -41/36 + 1/122*I, 4174/155 -
8/75*I, -89/186 + 164/113*I, 15/28 - 22/41*I, 50/117 - 22/123*I, -182/23 -
159/161*I, 73/97 - 13/147*I, 145/47 + 43/48*I, 247/98 - 5/64*I, 541/193 +
17/155*I, 284/193 + 8/29*I, -19/116 - 9/20*I, -101/23 - 70/69*I, -6/31 +
125/178*I, 359/117 + 30/49*I, 22/23 - 674/187*I, 373/199 - 3/8*I, 413/179 +
199/181*I, -3 + 118/161*I, -39/14 + 39/86*I, 5/79 + 85/167*I, -4115/138 +
27/31*I, -13/4 + 311/158*I, 13/159 - 4/5*I, 29/13 + 151/71*I, -89/51 + 7/5*I,
41/35 - 483/95*I, 11/90 - 137/188*I, 587/104 + 107/64*I, 659/97 + 31/88*I,
-137/162 - 575/122*I, -103/178 - 242/129*I, 169/121 - 25/48*I, -134/53 +
363/47*I, -112/45 + 2/51*I, 191/138 - 41/104*I, 2 - 1004/53*I, -5/7 +
281/101*I, 1601/200*I]}


\\ this one makes magma or me fail, seems to be magma
{ptsBug= [139/45 - 35/12*I, -113/92 + 248/95*I, 199/100 + 185/18*I, 50/33 -
213/74*I, 295/57 + 74/67*I, 60/47 + 99/100*I, -61/20 - 329/16*I, -27/98 +
1613/78*I, -101/89 + 2/3*I, -46/63 - 101/8*I, 455/81 + 299/100*I, 1/3 -
4/11*I, 208/45 + 179/11*I, -7/6 + 641/63*I, -617/58 - 227/69*I]}
{ptsBug2=[1/14 + 257/33*I, 77/32 + 106/77*I, 38/13 - 73/42*I, -17/95 -
233/27*I, 127/38 + 3/40*I, -1/5 - 7/58*I, -11/38 + 48/67*I, -73/54 -
31/38*I, 381/85 + 38/53*I, -103/39 + 51/32*I, -3/10 + 6/13*I, 208/63 +
70/43*I, -909/86 - 81/91*I, 16/27 + 5/22*I, -625/73 - 183/28*I]}
\\these remain bad with tree
{ptsBad=[195/47 - 79/9*I, 52/23 - 29/71*I, 16/17 - 1060/59*I, -1/100 + 19/66*I,
43/20 -17/39*I]}
{ptsBad2=[121/166 - 22/49*I, -10403/296 + 556/227*I, 97/429 + 9/128*I, 85/377 + 7/376*I, 973/379 + 5/42*I];}
{ptsBad3=[1/7 - 99/97*I, 1/4 + 41/25*I, 799/100 + 562/81*I, 79/98 - 117/71*I, -517/67 + 11/21*I, -27/11 + 1/53*I, -49/20 + 1/33*I]}
{ptsBad4=[7/12 + 3688/129*I, 133/180 - 265/274*I, 113/166 - 1/19*I, 123/158 +
71/204*I, -70/277 + 152/97*I, -21/89 + 457/293*I, -162/299 - 3/16*I]}

ptsRandom(g=2) = {
  vector(2*g+1,k,
   tan(3.1*(random(1.)-.5))+I*tan(3.1*(random(1.)-.5))
  );
}

/* *************************************************************************
 * speed and quality tests
 * ************************************************************************* */

fail(descr,args) = return(print(descr,args)); \\ and return 0, thanks to gp !

/* multi square root */
speed_sqrt(expand=0,n=3*10^4,size=0) = {
  for(i=1,n,
  A=ptsRandom(if(size,size,random(6)+1));
  z = ptsRandom(0)[1];
  sqrt_affinereduction(A,z,expand)
  );
}
/* it is very surprising that limiting the number of square roots
   gains nothing until 100 digits, no matter the number of roots
 */

test_sqrt(n=10^6,size=0) = {
  for(i=1,n,
      A=ptsRandom(if(size,size,random(6)+1));
      z = ptsRandom(0)[1];
      if(real((ps=sqrt_affinereduction(A,z,1))/
          (sa=sqrt_affinereduction(A,z)))<=0,
        print("prod_sqrt/sqrt_aff=",ps/sa);
        print("BUG A=",A,"\n z=",z);
        return(0));
     );
   \\ test passed
   return(1);
}

/* spanning tree */

speed_spanning(n=10^3,g=2,output_bad=1) = {
  inf = 3;
  bad=[];
  for(k=1,n,
    A = ptsRandom(g);
    tree = max_spanning(A);
    if(tree[2]<inf,
      inf=tree[2];bad=bestappr(A,40);
      if(output_bad,printf("tau=%1.2f for\n   A = %s\n",inf,bad))
      );
    );
    printf("min tau = %1.3f\n",inf);
}
\\ 1.5ms for g=2
\\ 10ms for g=4
\\ 97ms for g=8
\\ most of time spent in asinh/atanh
\\ mintau = 0.2 < Pi/9 = 0.34
\\ so might be useful to split integral
/* to compare : the following is only
   10 times faster but output is
   easily 10^5 times worse
 */
tau_max_nco(A) = {
  my(tau,maxmod,diameter);
  tau = Pi/2;
  for(ii=1,#A-1,
      for(jj=1,#A,
        if(jj==ii||jj==ii+1,next());
        tau = min(tau,tau_3(A[ii],A[ii+1],A[jj]));
        \\ diameter FIXME
        );
     );
  return(0.95*tau);
}
speed_tau_chain(n=10^3,g=2) = {
  inf = 3;
  for(k=1,n,
    A = ptsRandom(g);
    inf = min(inf,tau_max_nco(A));
    );
   print("min tau =",inf);
}

/** graph tree and intersections **/



/* using gp */
plot_tree(A,title=1) = {
  my(X,Y,E,t);
  [E,t] = max_spanning(A);
  X = real(A); Y = imag(A);
  my(w,s);
  w=0; s = plothsizes();
  plotinit(w, s[1]-1, s[2]-1);
  my(a,ma); a = .1; ma = [1+a,-a;-a,1+a];
  [xm,xM] = [vecmin(X),vecmax(X)]*ma;
  [ym,yM] = [vecmin(Y),vecmax(Y)]*ma;
  plotscale(w, xm,xM,ym,yM);
  \\plotpointsize(w,10); \\ only for gnuplot...
  for(i=1,#E,
      my(a,b);
      [a,b]=E[i];
      \\plotpoints(w,X,Y);
      plotlines(w,[X[a],X[b]],[Y[a],Y[b]]);
      );
  plotmove(w,(xm+xM)/2,yM);
  if(title,plotstring(w,Strprintf("tau = %1.2f",t),9));
  plotclip(w);
  plotdraw([w,0,0]);
  t;
}

/* good and bad trees */
bad_trees(n=1000,d=5) = {
    vecsort(vector(n,k,
                A=vector(d,j,random(2.)-1+I*(random(2.)-1));
                concat(max_spanning(A),[A])),2);
}


/* plot using gnuplot */
writeedge(file,e,f) = {
    write1(file,Strprintf("%1.2f\t%1.2f\n",real(e),imag(e)));
    write1(file,Strprintf("%1.2f\t%1.2f\n",real(f),imag(f)));
    /* when plotting with vectors (plot "file" w vectors) */
    \\write1(file,Strprintf("%1.2f\t%1.2f\t%1.2f\t%1.2f\n",real(e),imag(e),real(f-e),imag(f-e)));
    write(file,"");
}
write_gnuplot_tree(C,labels=1,edgesfile="edges.data",labelfile="labels.data") = {
  my(int,a,b,xmin,xmax,ymin,ymax);
  my(A=C[iRoots]);
  system(Str("rm -f ",edgesfile));
  if(labels,
      system(Str("rm -f ",labelfile));
    );
  my(edges = C[iTree]);
  \\my(int = C[iIntersection]);
  xmin = xmax = 0; ymin = ymax = 0;
  for(i=1,#edges,
      edge = edges[i];
      a = A[edge[1]];
      b = A[edge[2]];
      writeedge(edgesfile,a,b);
      if(labels,
        xa = real(a); ya = imag(a);
        xb = real(b); yb = imag(b);
        xmin = min(xmin,min(xa,xb));
        xmax = max(xmax,min(xa,xb));
        ymin = min(ymin,min(ya,yb));
        ymax = max(ymax,min(ya,yb));
        /* show the edge */
        write1(labelfile,
          Strprintf("set arrow from %1.2f,%1.2f to %1.2f,%1.2f\n",xa,ya,xb,yb));
        /* and print its number */
        write1(labelfile,
          Strprintf("set label \"%i\" at %1.2f,%1.2f\n",i,(xa+xb)/2,(ya+yb)/2));
        );
     );
}

hcShowTree(C,edgesfile="edges.data",labelfile="labels.data") = {
  system("rm -f gnuplottree");
  write("gnuplottree",
      Strprintf("load \"%s\"\nplot \"%s\" w lp\n",labelfile,edgesfile));
  write_gnuplot_tree(C);
  system("gnuplot -persist gnuplottree");
}
  

/* graph of AbelJacobi Map */
AJgnuplot(C,xmin=-5,xmax=5,ymin=-5,ymax=5,n=50) = {
  my(x,y,z);
  system("rm -f AbJac.data");
  file = "AbJac.data";
  \\ just for integration init
  for(i=1,n,
    x = xmin+(xmax-xmin)*i/n;
    for(j=1,n,
    y = ymin+(ymax-ymin)*j/n;
    z=hcAbelJacobi(C,x+I*y)[1];
    write1(file,Strprintf("%1.3f\t%1.3f\t%1.3f\t%1.3f\n",x,y,real(z),imag(z)));
    );
    write1(file,"\n");
  );
  system("gnuplot -persist AbJac.gnuplot");
}
  

/* *************************************************************************
 * some tests based on wanted properties
 */

test_intersection_signs(n=100) = {
  default(debug,1);
  for(i=1,n,
      pts = ptsRandom(1+random(7));
      tree = max_spanning(pts);
      intersection_spanning(A,tree[1]);
     );
}

/* test if Omega is symetric to current prec */
Riemann_symmetry(Omega) = {
  my(Sym);
  Sym = abs(Omega - mattranspose(Omega));
  if(trace(Sym*mattranspose(Sym))>1e-20,
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
Riemann_positivity(Omega,n=5) = {
  my(IOmega=imag(Omega));
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
  my(Omega);
  /* test if the small period matrix is symetric */
  C = hcInit(pts);
  Omega = hcSmallperiods(C);
  if(!Riemann_symmetry(Omega),return(fail("symmetry : ",pts)));
  if(!Riemann_positivity(Omega),return(fail("positivity : ",pts)));
  return(1);
}


test_Riemannrelations_random(n=100) = {
  for(i=1,n,
    pts = ptsRandom(random(7)+1);
    test_Riemannrelations(pts) || return(0);
    );
  return(1);
}


/* *************************************************************************
 *  compare using magma
 */

/* find an integral change of basis between U and V, based on the first line by
 default (or specify line i)
 returns P such that U*P = V 
*/
integral_base_change(U,V,prec=10,line=1) = {
  my(g,U1,res,j,cstprec,tmp);
  if(type(U)=="t_MAT",
    g = matsize(U)[1];
    U = U[line,],
    g = #U/2
    );
  if(type(V)=="t_MAT",
    ncol = min(#V[line,],2*g);
    V = V[line,];
    V = vecextract(V,Str("1..",ncol))
  );
  res = matrix(#V,2*g);
  U1 = vector(2*g+1);
  for(j=1,2*g,U1[j] = U[j]);
  /* do a lindep with low precision */
  cstprec = 1.+0.e0*10^-prec;
   for(i=1,#V,
    U1[2*g+1] = V[i];
    tmp = lindep(U1*cstprec);
    if(tmp[2*g+1] == 1,
      \\OK
      for(j=1,2*g,res[j,i] = -tmp[j])
      ,if(tmp[2*g+1]==-1,
        \\ still OK
        for(j=1,2*g,res[j,i] = tmp[j])
        ,
        print("no integral base change found (",tmp[2*g+1],"!=1)");
        return(0);
        );
      );
    );
   return(res); \\ should be a matrix in Sl(2g,Z)
}

matJ(g)=matrix(2*g,2*g,i,j,if(i+g==j,1,if(i==j+g,-1,0)));

/* check if a matrix is symplectic : should then return
   the canonical symplectic matrix.
   In particular, if the input is the change of basis between a symplectic
   period matrix (the result of magma) and a matrix U, it returns
   the intersection matrix corresponding to U
*/
intersection_from_change(change) = {
  g = matsize(change)[1]/2;
  M = change^(-1);  
  Jg = matJ(g);
  res = matrix(2*g,2*g,i,j,M[i,]*Jg*(M[j,]~));
  res
}

 /* useful for reading magma results, output via
   PrintFile("tmp.magma",$1,"Magma":Overwrite:=true);
 */
read_magma(file="tmp.magma") = {
  extern(Str("sed 's/.*\\[.*|/[/;' ../../magma/",file," | tr -d '\n\\\\'"))
}

/* compute period matrix with magma */
str_magmacommand(A,Dcalc=60,time=0,BigorSmall="Big",file="tmp.magma") = {
if(time,timestr="print \"  Time magma : \",t;\n",timestr="");
Strprintf(
"C<I> := ComplexField(%i);
K<x> := PolynomialRing(C);
function periods(rts)
  /* f(x) := prod_{r in rts}(x-r) */
  g := (#rts-1) div 2;
  f := &*[x-a : a in rts];
  fx := Evaluate(f,x);
  A := AnalyticJacobian(fx);
  M := %sPeriodMatrix(A);
  return M;
end function;
A:=%Ps;
t := Cputime();
X:=periods(A);
t := Cputime(t);
PrintFile(\"%s\",X,\"Magma\":Overwrite:=true);
%squit;
",Dcalc,BigorSmall,A,file,timestr);
}
magma_run(command) = {
  if(!system("which -s magma"),error("cannot find magma"));
  system("rm -f ../../magma/inputmagma");
  write("../../magma/inputmagma",command);
  \\ -b removes junk output, no need to > /dev/null
    system("cd ../../magma && magma -b inputmagma");
}
magma_Bigperiods(A,Dcalc=60,time=0) = {
  magma_run(str_magmacommand(A,Dcalc,time,"Big"));
  V = read_magma();
  g = floor((#A-1)/2);
  V = matrix(g,2*g,i,j,V[2*g*(i-1)+j]);
  return(V);
}
magma_Smallperiods(A,Dcalc=60,time=0) = {
  magma_run(str_magmacommand(A,Dcalc,time,"Small"));
  V = read_magma();
  g = floor((#A-1)/2);
  V = matrix(g,g,i,j,V[g*(i-1)+j]);
  return(V);
}

basechange_gp_magma(A,U=[]) = {
  V = magma_Bigperiods(A);
  if(!U,
    U = hcBigperiods(hcInit(A));
  );
  /* change of basis between me and magma */
  return(integral_base_change(U,V[1,]));
}

test_intersection_matrix(A,returndiff=1) = {
  C = hcInit(A);
  U = hcTreeperiods(C);
  my_int = hcIntersection(C);
  M = basechange_gp_magma(A,U);
    /* deduce intersection matrix of my tree basis, should be equal to Int_matrix*/
  magma_int = intersection_from_change(M);
  tmp = matsize(my_int);
  if(returndiff,
      return(matrix(tmp[1],tmp[2],i,j,Str("g",my_int[i,j],"/m",magma_int[i,j])))
      ,
      for(i=1,tmp[1],for(j=1,tmp[2],
          if(my_int[i,j]!=magma_int[i,j],
            print("## Erreur signe ",i,".",j," :
              g",my_int[i,j],"/m",magma_int[i,j]);
            return(0););
          ););
      return(1);
    );
}

/* test_intersection_matrix(pts=bestappr(ptsRandom(1),100)) */

test_intersection_matrix_random(n=200) = {
  i=0;pts=[];
  while(i++<=n&&test_intersection_matrix(pts=bestappr(ptsRandom(1+random(7)),100),0),);
  if(i<=n,print(pts);return(0));
  return(1);
}

/* test accuracy */
matrixnorm(X)={sqrt(trace(conj(mattranspose(X))*X))}
test_accuracy_magma(A) = {
  Dcalc=precision(1.);
  V = magma_Bigperiods(A,Dcalc,time=1);
  U = hcBigperiods(hcInit(A));
  P = integral_base_change(U,V[1,]);
  matrixnorm(V-U*P);
  \\V-U*mattranspose(P)
  }

/* test relations for magma implantation */
test_Riemannrelations_magma(pts) = {
  my(Omega);
  /* test if the small period matrix is symetric */
  Omega = magma_Smallperiods(pts,,1);
  if(!Riemann_symmetry(Omega),return(fail("symmetry : ",pts)));
  if(!Riemann_positivity(Omega),return(fail("positivity : ",pts)));
  return(1);
}
test_Riemannrelations_magma_random(n=100) = {
  for(i=1,n,
    pts = bestappr(ptsRandom(random(7)+1),100);
    test_Riemannrelations_magma(pts) || return(0);
    );
  return(1);
}

/* *************************************************************************
 * compare with agm and Richelot 
 */

/* Gauss's AGM
   I = int_a^b dx/sqrt((x-a)(x-b)(x-c))
     = Pi/agm(sqrt(c-a),sqrt(c-b))
*/
gauss_agm(P) = {
  a = P[1]; b = P[2] ; c = P[3];
  Pi/agm(sqrt(c-a),sqrt(c-b))
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
  epsilon=10^-(default(realprecision));
  while(1,
    A=sqrt((c-b)*(c-bp)*(cp-b)*(cp-bp));
    B=sqrt((c-a)*(c-ap)*(cp-a)*(cp-ap));
    C=sqrt((b-a)*(b-ap)*(bp-a)*(bp-ap));
    dba = b+bp-a-ap;
    dcb = c+cp-b-bp;
    dca = c+cp-a-ap;
    a1 = (c*cp-a*ap-B)/dca;
    a1p= (b*bp-a*ap-C)/dba;
    b1 = (b*bp-a*ap+C)/dba;
    b1p= (c*cp-b*bp-A)/dcb;
    c1 = (c*cp-b*bp+A)/dcb;
    c1p= (c*cp-a*ap+B)/dca;
    D1 = 4*D*(a*ap*dcb-b*bp*dca+c*cp*dba)/(dba*dcb*dca);
    a=a1; ap=a1p; b=b1; bp=b1p; c=c1; cp=c1p; D=D1;
    /* test */
    if(abs(a-ap)<epsilon && abs(b-bp)<epsilon && abs(c-cp)<epsilon,
			return(Pi*sqrt(D)*(S[1]*a+S[2])/((a-b)*(a-c)))
    )
  );
}


/* period([a,b,c]) ~ agm()
   period([a,b,c,d,e,f])~richelot([a,b,c,d,e,f],[0,1])
*/
compare(P) = {
  local(g,p,i);
  if(#P==3,
    gettime();
    g = gauss_agm(P);
    print("temps agm : ",gettime(),"ms.");
    p = period_old(P);
    print("temps integration : ",gettime(),"ms.");
    print("erreur : ",g-p);
    return(g)
    );
  if(#P==6,
    gettime();
    g = richelot(P,[0,1]);
    print("temps richelot : ",gettime(),"ms.");
    p = period_old(P);
    print("temps integration : ",gettime(),"ms.");
    print("erreur : ",g-I*p);
    i = period_intnum(P);
    print("temps pari intnum : ",gettime(),"ms.");
    print("erreur : ",g-I*i);
    return(g)
    );
  return(period_old(P));
}

/* *************************************************************************
 * automatic set of tests
 */
do_tests(nrels=30,nmagma=10,nsqrt=1000) = {
    test_sqrt(nsqrt)
 && !print("# sqrt passed")
 && test_Riemannrelations(ptsA)
 && test_Riemannrelations(ptsHexa)
 && test_Riemannrelations_random(nrels)
 && !print("# Riemann relations passed")
 /* is the intersection matrix correct ? */
 && test_intersection_matrix(ptsA,0)
 && test_intersection_matrix(ptsHexa,0)
 && test_intersection_matrix_random(nmagma)
 && !print("# intersection passed")
}
  
doplot(m=1) = {
  C = hcInit(ptsEx2);
  \\ploth(x=-0.5,1.5,z=hcAbelJacobi(C,(5/2+3/2*I)*x);[z[1],1-z[1]],,50)
  ploth(x=-0.5,1.5,z=hcAbelJacobi(C,(5/2+3/2*I)*x);z[1],,50)
  \\ploth(x=-0.5,1.5,z=hcAbelJacobi(C,(5/2+3/2*I)*x);[z[1],z[2],z[3]],,50);
}
domoebius(m=[1,2;2,3]) = {
  C1 = hcInit(ptsEx2);
  Omega1 = hcBigperiods(C1);
  C2 = hcInit(homog(m,ptsEx2));
  Omega2 = hcBigperiods(C2);
  alpha = moebius_matrix(m,C1[iGenus]);
}
/* We have lost signs of loops between intersection matrices I1 and I2 : chose
 * signs such that P^t*I1*P=I2 */
sign_changes(I1,I2) = {
  my(g = matsize(I1)[1]/2);
  my(diag=vector(2*g,i,0));
  my(Ip=I1);
  my(P=matid(2*g));
  diag[1]=1;
  i=1;j=2;
  \\while(Ip!=I2,
    \\ is there a problem at [i,j] ?
  \\   );
  for(i=1,2*g,
    for(j=i+1,2*g,
      if(Ip[i,j]!=I2[i,j],
          \\print("diff on coeff[",i,",",j,"]");
        if(diag[i]&&diag[j],return("Impossible"));
        if(diag[i],
                   \\print("i=",i," fixed, change j=",j);
                   diag[j]=-1;
                   P[j,j]=-1;
                   Ip[,j] *= -1; Ip[j,] *= -1
        ,
                   \\print("change i=",i);
                   diag[i]=-1;
                   P[i,i]=-1;
                   Ip[,i] *= -1; Ip[i,] *= -1;
         );
         );
     );
     \\if(!diag[i],diag[i]=1);
     );
   return(P);
}
