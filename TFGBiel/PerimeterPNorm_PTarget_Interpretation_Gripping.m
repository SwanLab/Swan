close all;
clear;

load('Gripping.mat');
clc;

mesh           = createMesh(mesh);
designVariable = createDesignVariable(mesh,fun);
clear('fun');
designVariable.fun.plot();

% Continue here






function mesh = createMesh(m)
    s.connec = m.connec;
    s.coord  = m.coord;
    mesh     = Mesh.create(s);
end


function designVariable = createDesignVariable(mesh,f)
    sF.fValues     = f.fValues;
    sF.mesh        = mesh;
    sF.order       = 'P1';
    s.fun          = LagrangianFunction(sF);
    s.mesh         = mesh;
    s.type         = 'Density';
    s.plotting     = false;
    dens           = DesignVariable.create(s);
    designVariable = dens;
end