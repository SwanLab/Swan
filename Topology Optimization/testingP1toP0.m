%% Create sample FEM results
clear; close all;

file = 'test2d_triangle';
a.fileName = file;
s = FemDataContainer(a);
fem = FEM.create(s);
fem.solve();

uCol   = fem.variables.d_u;
u      = reshape(uCol,[s.mesh.ndim,s.mesh.nnodes])';

%% P1 to P0 v2
cc.mesh   = s.mesh;
cc.connec = s.mesh.connec;
cc.nelem  = size(s.mesh.connec,1);
cc.nnode  = size(s.mesh.connec,2);
cc.npnod  = size(s.mesh.coord,1);
projector2 = Projector_P1toP0(cc);

% Create a FeFunction
z.connec = s.mesh.connec;
z.type   = s.mesh.type;
z.fNodes = u;
uFeFun = P1Function(z);

u_P0 = projector2.project(uFeFun);

% Plot using new stuff
figure()
[m, f] = createDiscontP0(s.mesh, u_P0(:,1));
plotDiscontP0(m,f);

% Plot using old stuff
figure()
aa.connec = s.mesh.connec;
aa.type   = s.mesh.type;
aa.fNodes = u;
fefunDisp = P1Function(aa);
p0displac = fefunDisp.computeValueInCenterElement()';
plotDiscontP0(m,f);

% Plot using P0Function
pp.mesh   = s.mesh;
pp.fElem = u_P0(:,1);
P0Fun     = P0Function(pp);
P0Fun.plot();

function [meshDisc, fDisc] = createDiscontP0(mesh,f)
    meshDisc = mesh.createDiscontinousMesh();
    nnodeElem = meshDisc.nnodeElem;
    fRepeted = zeros(size(f,1),nnodeElem);
    for iNode = 1:nnodeElem

        fRepeted(:,iNode) = f;
    end
    fRepeted = transpose(fRepeted);
    fDisc = fRepeted;
end

function plotDiscontP0(mesh, f)
    trisurf(mesh.connec, mesh.coord(:,1), mesh.coord(:,2), f)
end 