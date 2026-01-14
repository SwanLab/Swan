function [xx, ww] = GaussQuad(n, a, b)
if nargin == 0
    n=20; 
    a = -1 ; 
    b=1 ;
end
%GAUSSQUAD Gaussian quadrature integration.
%
%   [X, W] = GAUSSQUAD(N, A, B) returns the base points X and weight factors W
%   for integrating a function from A to B using an N-point Gaussian quadrature
%   rule.  This integral is exact for a polynomial of degree 2*N-1 or lower.
%
%   [X, W] = GAUSSQUAD(N, C) assumes the interval is from 0 to C.
%GaussQuad
%   [X, W] = GAUSSQUAD(N) assumes the interval is from -1 to 1.
%
%   GAUSSQUAD(...) with no output arguments plots the base points X and the
%   weight factors W.

%   This program is based on an anonymous MATLAB program found on the MathWorks
%   FTP server, ftp://ftp.mathworks.com, extended and rewritten for speed.
%
%   According to the original program, concepts for finding the base points and
%   weight factors can be found on page 93 of "Methods of Numerical
%   Integration" by Philip Davis and Philip Rabinowitz.

%   Author:      Peter J. Acklam
%   Time-stamp:  2004-02-02 08:36:19 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   % Check number of input arguments.
%    narginchk(1, 3);
% 
%    % Assign default values to missing arguments.
%    switch nargin
%       case 1                    % GAUSSQUAD(N)
%          b = 1;
%          a = -1;
%       case 2                    % GAUSSQUAD(N, C)
%          b = a;
%          a = 0;
%    end

   u = 1 : n-1;
   u = u ./ sqrt(4*u.^2 - 1);

   % Same as A = diag(u, -1) + diag(u, 1), but faster (no addition).
   A = zeros(n, n);
   A( 2 : n+1 : n*(n-1) ) = u;
   A( n+1 : n+1 : n^2-1 ) = u;

   % Find the base points X and weight factors W for the interval [-1,1].
   [v, x] = eig(A);
   [x, k] = sort(diag(x));
   w = 2 * v(1,k)'.^2;

   % Linearly transform from [-1,1] to [a,b].
   x = (b - a) / 2 * x + (a + b) / 2;   % transform base points X
   w = (b - a) / 2 * w;                 % adjust weigths

   % If output arguments are given, return output and exit.
   if nargout
      xx = x;
      ww = w;
      return
   end

   % Plot base points and weight factors.
   userdata = 'Gaussian Quadrature Base Points and Weight Factors';
   fig = findobj(0, 'Type', 'figure', 'UserData', userdata);
   if isempty(fig)
      fig = figure('UserData', userdata);
   end
   figure(fig);
   h = plot(x, zeros(1, n), '.', x, w, '.');
   set(h, 'MarkerSize', 12);
   set(gca, 'XLim', [a b]);
   title(sprintf('Gaussian Quadrature with %d points', n));
   xlabel('Base points');
   ylabel('Weight factors');
