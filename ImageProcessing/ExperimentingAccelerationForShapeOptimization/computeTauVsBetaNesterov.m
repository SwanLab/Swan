function [tauGrid,betaGrid,iters,converged] = computeTauVsBetaNesterov()
    beta = 0:0.05:1;
    tau  = 10:10:200;
    [betaGrid,tauGrid] = meshgrid(beta,tau);
    iters = zeros(size(betaGrid,1),size(betaGrid,2));
    converged = iters;
    for i = 1:numel(tau)
        for j = 1:numel(beta)
            [iters(i,j), converged(i,j)] = computeShapeOptimization(betaGrid(i,j),tauGrid(i,j));
            fprintf('Beta: %.3f, Tau: %.3f, Iterations: %.0f\n',betaGrid(i,j),tauGrid(i,j),iters(i,j));
        end
    end
end

function [it,converged] = computeShapeOptimization(beta,tau)
    s.momentumParams.type = 'CONSTANT';
    s.momentumParams.value = beta;
    s.tau = tau;
    s.maxIter = 1000;
    s.TOL = 1e-3;
    solver = ShapeOptimizationSolver(s);
    solver.solve();
    it = numel(solver.incXvalues);
    converged = solver.incXvalues(end) < 1e-2;
    m = solver.designVariable.mesh;
    TO = triangulation(m.connec,m.coord(:,1),m.coord(:,2),solver.designVariable.value);
    f = figure();
    t = trimesh(TO);
    view(0,90)
    c = colorbar;
    colormap(flipud(gray))
    grid off
    filename = "./ImageProcessing/ExperimentingAccelerationForShapeOptimization/Cantilever_results/" + string(beta) + "_" + tau + ".png";
    saveas(f,filename);
    close(f);
end