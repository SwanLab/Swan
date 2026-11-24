function h = plotLatticeBW(dehomogenizer,invert,fileName)
if nargin < 2 || isempty(invert)
    invert = false;
end
if nargin < 3
    fileName = '';
end
mesh = dehomogenizer.getMesh();
psi  = dehomogenizer.getLevelSet().fValues;
colors = [1 1 1; 0 0 0];
if invert
    colors = flipud(colors);
end
h = figure('Color','w');
patch('Faces',mesh.connec,'Vertices',mesh.coord, ...
    'FaceVertexCData',psi,'FaceColor','interp','EdgeColor','none');
colormap(colors);
clim([-1e-9,1e-9]);
axis equal off
if ~isempty(fileName)
    print(h,fileName,'-dpng','-r300');
end
end