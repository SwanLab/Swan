
clear;
close all;
clc;

load('DVToStudyEpsilon.mat','d');
d.fun.plot();
a = gcf().findobj();
a(4).EdgeColor = 'none';

% Here create your filter of interest
filter = createAniFilter(d);

rhoEps = filter.compute(d.fun,2);
rhoEps.plot();
a = gcf().findobj();
a(4).EdgeColor = 'none';










% Functions definition
function aniFilter = createAniFilter(d)
    mesh = d.fun.mesh;
    h    = mesh.computeMeanCellSize();

    s.filterType   = 'PDE';
    s.mesh         = mesh;
    s.boundaryType = 'Neumann';
    s.metric       = 'Anisotropy';
    nu             = 45;  
    s.aniAlphaDeg  = 90;   
    epsilon        = h;
    s.CAnisotropic = [tand(nu), 0; 0, 1/tand(nu)];  
    s.trial        = LagrangianFunction.create(mesh, 1, 'P1');
    aniFilter      = Filter.create(s);
    aniFilter.updateEpsilon(epsilon);
end