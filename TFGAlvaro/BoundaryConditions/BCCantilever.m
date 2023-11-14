clear

% Cantilever BCs
%load('TFGAlvaro/BoundaryConditions/dataCantilever.mat','lnodes','pointload_complete')
run('TFGAlvaro/GiD/CantileverAlvaro.gid/CantileverAlvaro.m')

lnodes = [];
pointload_complete = [];
coor = gidcoord(:,2:end);

% Nodes constrained
