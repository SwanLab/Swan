params.designVariable = struct('thetavec', [1.0, 2.0, 3.0]);
params.x = double([1.0, 2.0, 3.0]);


output = callJuliaClass('Sh_Func_L2norm', 'computeStochasticCostAndGradient', params);
disp(output.j)   % cost
disp(output.dj)  % gradient