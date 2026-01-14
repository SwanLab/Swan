function hh = quadplot(quad,varargin)
%TRIPLOT Plots a 2D triangulation
%   QUADPLOT(QUAD,X,Y) displays the quadrilaterals defined in the
%   M-by-4 matrix QUAD.  A row of QUAD contains indices into X,Y that
%   define a single quadrilateal. The default line color is blue.
%
%   QUADPLOT(...,COLOR) uses the string COLOR as the line color.
%
%   H = QUADPLOT(...) returns a vector of handles to the displayed 
%   quadrilaterals
%
%   QUADPLOT(...,'param','value','param','value'...) allows additional
%   line param/value pairs to be used when creating the plot.
%
%   See also TRISURF, TRIMESH, DELAUNAY, TriRep, DelaunayTri.
%
%   Script code based on copyrighted code from mathworks for TRIPLOT.
%   Allan P. Engsig-Karup, apek@imm.dtu.dk.

error(nargchk(1,inf,nargin,'struct'));

start = 1;

x = varargin{1};
y = varargin{2};
quads = quad;
if (nargin == 3) || (mod(nargin-3,2) == 0)
    c = 'blue';
    start = 3;
else
    c = varargin{3};
    start = 4;
end
  
d = quads(:,[1 2 3 4 1])';
h = plot(x(d), y(d),c,varargin{start:end});
if nargout == 1, hh = h; end
