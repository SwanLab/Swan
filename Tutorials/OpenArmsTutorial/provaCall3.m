% Prepare input parameters for constructor
params.designVariable = struct('thetavec', [1.0, 2.0, 3.0]);

% Call the constructor method via the Julia bridge
output = JuliaShFuncL2norm(params);

% Display result
%disp(output.status);       % should be "created"
%disp(output.thetavec);     % should be [1.0 2.0 3.0]