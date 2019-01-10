install("hcinit","Gp","hcinit","./hyperellperiods.so");
install("hyperellperiods","GD0,L,p","hyperellperiods","./hyperellperiods.so");
C1 = hyperellperiods(x^3-x,1)
C2 = hyperellperiods(x^3-x,1)
/* periods of elliptic curves */
T(l)=hyperellperiods(x*(x-1)*(x-l),1);
T2(l)=my([w1,w2]=ellperiods(ellinit([0,-l-1,0,l,0])));w1/w2;[w1,w2];
L(t)={
  my(N=ceil(sqrt(precision(1.)*log(10)/Pi*imag(t)))+3);
  my(q=exp(I*Pi*t));
  (sum(n=-N,N,q^(n+1/2)^2)/sum(n=-N,N,q^n^2))^-4;
}
