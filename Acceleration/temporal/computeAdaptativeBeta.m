function [tau,iters] = computeAdaptativeBeta()
    tau = 50:10:200;
    iters = zeros(1,numel(tau));
    for i = 1:numel(tau)
        [it,converged] = computeShapeOptimization(tau(i));
        fprintf('Tau: %.3f, Iterations: %.0f\n',tau(i),it);
        iters(i) = it;
    end
    save adaptative_beta_bridge1 iters tau
end

function [it,converged] = computeShapeOptimization(tau)
    s.momentumParams.type = 'CONVEX';
    s.tau = tau;
    s.maxIter = 300;
    s.TOL = 1e-3;
    solver = ShapeOptimizationSolver(s);
    solver.solve();
    it = numel(solver.incJV);
    converged = solver.incJV(end) < 1e-3;
    m = solver.designVariable.mesh;
    TO = triangulation(m.connec,m.coord(:,1),m.coord(:,2),solver.designVariable.value);
    f = figure();
    t = trimesh(TO);
    view(0,90)
    c = colorbar;
    colormap(flipud(gray))
    grid off
    filename = "./ImageProcessing/ExperimentingAccelerationForShapeOptimization/Bridge_1_results/" + "ADAPTATIVE" + "_" + tau + ".png";
    saveas(f,filename);
    close(f);
end