function h = plotDesignRhoBW(designFile,fileName,invert)
if nargin < 2
    fileName = '';
end
if nargin < 3 || isempty(invert)
    invert = false;
end
d   = load(designFile);
rho = mean(d.rho(d.connec),2);
colors = flipud(gray);
if invert
    colors = gray;
end
h = figure('Color','w');
patch('Faces',d.connec,'Vertices',d.coord,'FaceVertexCData',rho, ...
    'FaceColor','flat','EdgeColor','none');
colormap(colors);
clim([0,1]);
axis equal off
if ~isempty(fileName)
    print(h,fileName,'-dpng','-r300');
end
end