function [p,S,mu] = polyfix(x,y,n,xfix,yfix,xder,dydx)
%POLYFIX Fit polynomial p to data, but specify value at specific points
%
%   p = polyfix(x,y,n,xfix,yfix) 
% finds the coefficients of the polynomial of degree n that fits the data 
% in a least-squares sense, with the constraint that polyval(p,xfix) = yfix
%
%   p = polyfix(x,y,n,xfix,yfix,xder,dydx)
% uses the additional constraint that the derivative at xder = dydx
%   
%   [p,S] = polyfix(...) or [p,S,mu] is also possible
%   See the documentation for polyfit for details
%
% The polynomial order n must be high enough to match all specified values
% NOTE: For the lowest order allowed, p will fit the constraints, but 
%       may disregard x and y.
%
% Example 1:
% x = linspace(0,2,100)';y = sin(pi*x)+ 0.1*randn(100,1);
% p = polyfix(x,y,3,[0,2],[0,0]);plot(x,y,'.',x,polyval(p,x));
% 
% Example 2:
% x = linspace(0,1,100)';y = sin(x*pi/2) + 0.1*randn(100,1);
% p = polyfix(x,y,4,[],[],[0 1],[1 0]);plot(x,y,'.',x,polyval(p,x))
% See also: polyfit, polyval
% Are Mjaavatten, Telemark University College, Norway, November 2015
% Revision history
% 2015-11-28: Version 1.0
% 2015-12-07: Version 1.1:
%             Added option for specifying derivatives
%             The output is now a row vector, for consistency with polyfit
%             Added test for polynomial degree
% 2017-08-23: Version 1.2
%             Fixed trivial errors in the help section
% 2019-10-30: Version 1.3
%             Added outputs S and mu, in line with polyfit
  %% Make sure all input arrays are column vectors of compatible sizes:
  x = x(:);
  y = y(:);
  % Center and scale if the user specifies mu as an output
  if nargout > 2
     mu = [mean(x); std(x)];
     x = (x - mu(1))/mu(2);
     xfix = (xfix - mu(1))/mu(2);
  end
  
  nfit = length(x);
  if ~(length(y)== nfit)
      error('x and y must have the same size');
  end
  xfix = xfix(:);
  yfix = yfix(:);
  nfix = length(xfix);
  if ~(length(yfix)== nfix)
      error('xfit and yfit must have the same size');
  end
  if nargin > 5 % Derivatives specified
      xder = xder(:);
      dydx = dydx(:);
      if nargout > 2
        xder = (xder-mu(1))/mu(2);
        dydx = dydx*mu(2);
      end     
  else
      xder = [];
      dydx = [];
  end
  nder = length(xder); 
  if ~(length(dydx) == nder)
      error('xder and dydx must have the same size');
  end
  nspec = nfix + nder;
  specval = [yfix;dydx];
  %% First find A and pc such that A*pc = specval
  A = zeros(nspec,n+1); 
  % Specified y values
  for i = 1:n+1
      A(1:nfix,i) = ones(nfix,1).*xfix.^(n+1-i);
  end
  % Specified values of dydx
  if nder > 0
      for i = 1:n
          A(nfix +(1:nder),i) = (n-i+1)*ones(nder,1).*xder.^(n-i);
      end
  end
  if nfix > 0
      lastcol = n+1;
      nmin = nspec - 1;
  else
      lastcol = n;   % If only derivatives, p(n+1) is arbitrary
      nmin = nspec;
  end
  if n < nmin
      error('Polynomial degree too low. Cannot match all constraints');
  end    
  %% Find the unique polynomial of degree nmin that fits the constraints. 
  firstcol = n-nmin+1;   % A(:,firstcol_lastcol) detrmines p0   
  pc0 = A(:,firstcol:lastcol)\specval;  % Satifies A*pc = specval
  % Now extend to degree n and pad with zeros:
  pc = zeros(n+1,1);
  pc(firstcol:lastcol) = pc0;    % Satisfies A*pcfull = yfix
  % Construct Vandermonde matrix.
  V(:,n+1) = ones(length(x),1,class(x));
  for j = n:-1:1
     V(:,j) = x.*V(:,j+1);
  end
  % Subtract constraints polynomial values from y. 
  yfit = y-polyval(pc,x);
  %% We now find the p0 that minimises (V*p0-yfit)'*(V*p0-yfit)
  %  given that A*p0 = 0
  B = null(A);    % For any (n+1-nspc by 1) vector z, A*B*z = 0     
  z = V*B\yfit;   % Least squares solution of V*B*z = yfit
  p0 = B*z;       % Satisfies A*p0 = 0;
  p = p0'+pc';    % Satisfies A*p = b; 
  %% Add error etimation struct in the same way as polyfit
  if nargout > 1
    [~,R] = qr(V,0);
    r = y - V*p';
    % S is a structure containing three elements: the triangular factor from a
    % QR decomposition of the Vandermonde matrix, the degrees of freedom and
    % the norm of the residuals.
    S.R = R;
    S.df = max(0,length(y) - (n+1));
    S.normr = norm(r);
  end
end      
