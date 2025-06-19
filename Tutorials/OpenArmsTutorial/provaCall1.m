params.designVariable = struct('thetavec', [1.0, 2.0, 3.0]);
params.x = [1.0, 2.0, 3.0];

output = callJuliaClass('Sh_Func_L2norm', 'computeFunctionAndGradient', params);
disp(output.j)   % cost
disp(output.dj)  % gradient