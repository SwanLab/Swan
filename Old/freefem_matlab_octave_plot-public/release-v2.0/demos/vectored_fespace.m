% Author: Chloros2 <chloros2@gmx.de>
% Created: 2019-02-27
%
% Copyright (C) 2018 Chloros2 <chloros2@gmx.de>
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see
% <https://www.gnu.org/licenses/>.
%

clear all;

addpath('ffmatlib');

[p,b,t,nv,nbe,nt] = ffreadmesh('vectored_mesh.msh');
%Vector valued space {'P2','P2','P1'}
[vhIn] = ffreaddata('vectored_vh.txt');
%Solutions u,v,w - they are all the same. Need only one component.
[uIn] = ffreaddata('vectored_data.txt');
%[uIn,vIn,wIn] = ffreaddata('vectored_data3.txt');

%Extract first component (PDE solution u and subspace for P2)
[vhP2,u] = ffvectorget({'P2','P2','P1'}, vhIn, uIn, 1);
%Extract second component - same subspace as for component #1
[vhP2,v] = ffvectorget({'P2','P2','P1'}, vhIn, uIn, 2);
%Extract third component
[vhP1,w] = ffvectorget({'P2','P2','P1'}, vhIn, uIn, 3);

ffpdeplot(p,b,t,'VhSeq',vhP2,'XYData',u,'Mesh','on','MColor','b','Boundary','on','ColorMap',jet);
ylabel('y');
xlabel('x');
figure;
ffpdeplot(p,b,t,'VhSeq',vhP2,'XYData',v,'Mesh','on','MColor','b','Boundary','on','ColorMap',jet);
ylabel('y');
xlabel('x');
figure;
ffpdeplot(p,b,t,'VhSeq',vhP2,'FlowData',[u,v],'FGridParam',[24, 24],'Mesh','off','Boundary','on','ColorBar','off');
ylabel('y');
xlabel('x');
axis tight;
figure;
ffpdeplot(p,b,t,'VhSeq',vhP1,'XYData',w,'Mesh','on','MColor','b','Boundary','on','ColorMap',jet);
ylabel('y');
xlabel('x');
