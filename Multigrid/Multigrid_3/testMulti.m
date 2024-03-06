%pv are coordinate nodes (x,y)
%hmax refers to the maximal length allowed during meshing element generation by pmesh.m
%vUp and vDown refers to the prescribed number of restrictions steps to coarsen our grid, and vice versa 
%with vDown i.e. interpolations back into a finer grid
%iref is a prescribed number of iterations

%mginit initializes the multigrid solver and returns a structure array named data, containing all problem data
%mgsolve performs actual multigrid iterations via recursion calls on vcycle.m function and returns residual history and 
%estimated solution u.
clc;clear;close all;

tolerance = 1e-10;
vDown = 2;
vUp = 2;
hmax = 0.5;

o = MultigridTesting3;
dataTotal = o.getdata;
bcTotal = o.getBC;
meshTotal = o.getMesh;
% [u, res] = mgsolve(dataTotal, vDown, vUp, tolerance, bcTotal, meshTotal);
% semilogy(res), hold on

%pv = [0,0; 2,0; 1.5,1; .5,1; 0,0];
for iref = 5:-1:1
    tic
    % data = mginit(pv, hmax, iref); toc
    for i = 1:5-iref+1
        data(i) = dataTotal(iref+i-1);
        bc(i) = bcTotal(iref+i-1);
        mesh(i) = meshTotal(iref+i-1);
    end
    [u, res] = mgsolve(data, vDown, vUp, tolerance, bc, mesh);
    uTotal(:,5-iref+1) = u;
    a(5-iref+1) = toc;
    semilogy(res), hold on
    title('iter vs res')
    xlabel('iter')
    ylabel('res')
    legend('1 Mesh', '2 Meshes', '3 Meshes', '4 Meshes', '5 Meshes')
end
hold off
figure
plot(a)