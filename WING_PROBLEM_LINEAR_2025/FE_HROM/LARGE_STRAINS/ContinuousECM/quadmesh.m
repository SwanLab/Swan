function hh = quadmesh(quad,x,y,z,varargin)
%QUADMESH Quadrilateral mesh plot.
%   QUADMESH(QUAD,X,Y,Z,C) displays the quadrilaterals defined in the M-by-4
%   face matrix QUAD as a mesh.  A row of QUAD contains indexes into
%   the X,Y, and Z vertex vectors to define a single quadrilateral face.
%   The edge color is defined by the vector C.
%
%   QUADMESH(QUAD,X,Y,Z) uses C = Z, so color is proportional to surface
%   height.
%
%   QUADMESH(TRI,X,Y) displays the quadrilaterals in a 2-d plot.
%
%   H = QUADMESH(...) returns a handle to the displayed quadrilaterals.
%
%   QUADMESH(...,'param','value','param','value'...) allows additional
%   patch param/value pairs to be used when creating the patch object. 
%
%   See also PATCH.
%
% Script code based on copyrighted code from mathworks for TRIMESH.
% Allan P. Engsig-Karup, apek@mek.dtu.dk.

ax = axescheck(varargin{:});
ax = newplot(ax);

if nargin == 3 || (nargin > 4 && ischar(z))
  d = tri(:,[1 2 3 4 1])';
  if nargin == 3
    h = plot(ax, x(d), y(d));
  else
    h = plot(ax, x(d), y(d),z,varargin{1},varargin{2:end});
  end
  if nargout == 1, hh = h; end
  return;
end

start = 1;
if nargin>4 && rem(nargin-4,2)==1
  c = varargin{1};
  start = 2;
elseif nargin<3
  error(id('NotEnoughInputs'),'Not enough input arguments');
else
  c = z;
end

if ischar(get(ax,'color')),
  fc = get(gcf,'Color');
else
  fc = get(ax,'color');
end

h = patch('faces',quad,'vertices',[x(:) y(:) z(:)],'facevertexcdata',c(:),...
	  'facecolor',fc,'edgecolor',get(ax,'defaultsurfacefacecolor'),...
	  'facelighting', 'none', 'edgelighting', 'flat',...
      'parent',ax,...
	  varargin{start:end});
if ~ishold(ax), view(ax,3), grid(ax,'on'), end
if nargout == 1, hh = h; end

function str = id(str)
str = ['MATLAB:quadmesh:' str];

return

% % Example: plot quadmesh
% close all, clear all, clc
% 
% x = [0 1 2 2 1 1 0 0];
% y = [0 0 0 1 1 2 2 1];
% for i = 1 : length(x)
%     text(x(i)+0.1,y(i)+0.1,sprintf('%d',i))
%     hold on
% end
% quad = [1 2 5 8;
%     2 3 4 5;
%     8 5 6 7];
% 
% plot(x,y,'k.-')
% axis off
% 
% z = sin(pi*x).*cos(pi*y);
% c = z;
% figure
% 
% ax=newplot;
% fc = get(gcf,'Color');
% h = patch('faces',quad,'vertices',[x(:) y(:) z(:)],'facevertexcdata',c(:),...
%     'facecolor',fc,'edgecolor',get(ax,'defaultsurfacefacecolor'),...
%     'facelighting', 'none', 'edgelighting', 'flat',...
%     'parent',ax);
