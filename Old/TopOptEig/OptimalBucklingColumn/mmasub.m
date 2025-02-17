%    This is the file mmasub.m
%
function [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2, ...
f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a,c,d);
%
%    Written in May 1999 by
%    Krister Svanberg <krille@math.kth.se>
%    Department of Mathematics
%    SE-10044 Stockholm, Sweden.
%
%    This function mmasub performs one MMA-iteration, aimed at
%    solving the nonlinear programming problem:
%         
%      Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )
%    subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m
%                xmin_j <= x_j <= xmax_j,    j = 1,...,n
%                z >= 0,   y_i >= 0,         i = 1,...,m
%*** INPUT:
%
%   m    = The number of general constraints.
%   n    = The number of variables x_j.
%  iter  = Current iteration number ( =1 the first time mmasub is called).
%  xval  = Column vector with the current values of the variables x_j.
%  xmin  = Column vector with the lower bounds for the variables x_j.
%  xmax  = Column vector with the upper bounds for the variables x_j.
%  xold1 = xval, one iteration ago (provided that iter>1).
%  xold2 = xval, two iterations ago (provided that iter>2).
%  f0val = The value of the objective function f_0 at xval.
%  df0dx = Column vector with the derivatives of the objective function
%          f_0 with respect to the variables x_j, calculated at xval.
% df0dx2 = Column vector with the non-mixed second derivatives of the
%          objective function f_0 with respect to the variables x_j,
%          calculated at xval. df0dx2(j) = the second derivative
%          of f_0 with respect to x_j (twice).
%          Important note: If second derivatives are not available,
%          simply let df0dx2 = 0*df0dx.
%  fval  = Column vector with the values of the constraint functions f_i,
%          calculated at xval.
%  dfdx  = (m x n)-matrix with the derivatives of the constraint functions
%          f_i with respect to the variables x_j, calculated at xval.
%          dfdx(i,j) = the derivative of f_i with respect to x_j.
%  dfdx2 = (m x n)-matrix with the non-mixed second derivatives of the
%          constraint functions f_i with respect to the variables x_j,
%          calculated at xval. dfdx2(i,j) = the second derivative
%          of f_i with respect to x_j (twice).
%          Important note: If second derivatives are not available,
%          simply let dfdx2 = 0*dfdx.
%  low   = Column vector with the lower asymptotes from the previous
%          iteration (provided that iter>1).
%  upp   = Column vector with the upper asymptotes from the previous
%          iteration (provided that iter>1).
%  a0    = The constants a_0 in the term a_0*z.
%  a     = Column vector with the constants a_i in the terms a_i*z.
%  c     = Column vector with the constants c_i in the terms c_i*y_i.
%  d     = Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
%     
%*** OUTPUT:
%
%  xmma  = Column vector with the optimal values of the variables x_j
%          in the current MMA subproblem.
%  ymma  = Column vector with the optimal values of the variables y_i
%          in the current MMA subproblem.
%  zmma  = Scalar with the optimal value of the variable z
%          in the current MMA subproblem.
%  lam   = Lagrange multipliers for the m general MMA constraints.
%  xsi   = Lagrange multipliers for the n constraints alfa_j - x_j <= 0.
%  eta   = Lagrange multipliers for the n constraints x_j - beta_j <= 0.
%   mu   = Lagrange multipliers for the m constraints -y_i <= 0.
%  zet   = Lagrange multiplier for the single constraint -z <= 0.
%   s    = Slack variables for the m general MMA constraints.
%  low   = Column vector with the lower asymptotes, calculated and used
%          in the current MMA subproblem.
%  upp   = Column vector with the upper asymptotes, calculated and used
%          in the current MMA subproblem.
%
epsimin = sqrt(m+n)*10^(-9);
feps = 0.000001;
asyinit = 0.2;
asyincr = 1.1;
asydecr = 0.65;
albefa = 0.1;
een = ones(n,1);
zeron = zeros(n,1);

% Calculation of the asymptotes low and upp :
if iter < 2.5
  low = xval - asyinit*(xmax-xmin);
  upp = xval + asyinit*(xmax-xmin);
else
  zzz = (xval-xold1).*(xold1-xold2);
  factor = een;
  factor(find(zzz > 0)) = asyincr;
  factor(find(zzz < 0)) = asydecr;
  low = xval - factor.*(xold1 - low);
  upp = xval + factor.*(upp - xold1);
end

% Calculation of the bounds alfa and beta :
zzz = low + albefa*(xval-low);
alfa = max(zzz,xmin);
zzz = upp - albefa*(upp-xval);
beta = min(zzz,xmax);

% Calculations of p0, q0, P, Q and b.

ux1 = upp-xval;
ux2 = ux1.*ux1;
ux3 = ux2.*ux1;
xl1 = xval-low;
xl2 = xl1.*xl1;
xl3 = xl2.*xl1;
ul1 = upp-low;
ulinv1 = een./ul1;
uxinv1 = een./ux1;
xlinv1 = een./xl1;
uxinv3 = een./ux3;
xlinv3 = een./xl3;
diap = (ux3.*xl1)./(2*ul1);
diaq = (ux1.*xl3)./(2*ul1);
p0 = zeron;
p0(find(df0dx > 0)) = df0dx(find(df0dx > 0));
p0 = p0 + 0.001*abs(df0dx) + feps*ulinv1;
p0 = p0.*ux2;
q0 = zeron;
q0(find(df0dx < 0)) = -df0dx(find(df0dx < 0));
q0 = q0 + 0.001*abs(df0dx) + feps*ulinv1;
q0 = q0.*xl2;
dg0dx2 = 2*(p0./ux3 + q0./xl3);
del0 = df0dx2 - dg0dx2;
delpos0 = zeron;
delpos0(find(del0 > 0)) = del0(find(del0 > 0));
p0 = p0 + delpos0.*diap;
q0 = q0 + delpos0.*diaq;
P = zeros(m,n);
P(find(dfdx > 0)) = dfdx(find(dfdx > 0));
P = P * diag(ux2);
Q = zeros(m,n);
Q(find(dfdx < 0)) = -dfdx(find(dfdx < 0));
Q = Q * diag(xl2);
dgdx2 = 2*(P*diag(uxinv3) + Q*diag(xlinv3));
del = dfdx2 - dgdx2;
delpos = zeros(m,n);
delpos(find(del > 0)) = del(find(del > 0));
P = P + delpos*diag(diap);
Q = Q + delpos*diag(diaq);
b = P*uxinv1 + Q*xlinv1 - fval ;

%%% Solving the subproblem by a primal-dual Newton method
[xmma,ymma,zmma,lam,xsi,eta,mu,zet,s] = ...
subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d);




