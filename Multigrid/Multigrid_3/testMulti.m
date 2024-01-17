%pv are coordinate nodes (x,y)
%hmax refers to the maximal length allowed during meshing element generation by pmesh.m
%vUp and vDown refers to the prescribed number of restrictions steps to coarsen our grid, and vice versa 
%with vDown i.e. interpolations back into a finer grid
%iref is a prescribed number of iterations

%mginit initializes the multigrid solver and returns a structure array named data, containing all problem data
%mgsolve performs actual multigrid iterations via recursion calls on vcycle.m function and returns residual history and 
%estimated solution u.

tolerance = 1e-10;
vDown = 2;
vUp = 2;
hmax = 0.5;

pv = [0,0; 2,0; 1.5,1; .5,1; 0,0];
for iref = 1:3
    tic
    data = mginit(pv, hmax, iref); toc
    [u, res] = mgsolve(data, vDown, vUp, tolerance);
    semilogy(res), hold on
end
hold off