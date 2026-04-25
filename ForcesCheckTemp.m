%% Load data
S1 = load("RHSListTuron.mat");   % adapt name if needed
S2 = load("RHSListNoTuron.mat"); 

Sanalisys = S2;

FList = Sanalisys.RHSList;

nIter = size(FList,2);
ndof  = size(FList,1);

fprintf("==== RHS ANALYSIS START ====\n\n");

%% 1. Norm evolution
fprintf("---- RHS Norm ----\n");

Fnorm = zeros(nIter,1);

for i = 1:nIter
    F = FList(:,i);

    Fnorm(i) = norm(F);

    fprintf("Iter %2d → ||F|| = %.4e\n", i, Fnorm(i));
end

%% 2. Increment of RHS (key for instabilities)
fprintf("\n---- RHS Increment ----\n");

for i = 2:nIter
    dF = norm(FList(:,i) - FList(:,i-1));
    fprintf("Iter %2d → dF = %.4e\n", i, dF);
end

%% 3. Component-wise activity (which DOFs are driving it)
fprintf("\n---- Max component per iteration ----\n");

maxComp = zeros(nIter,1);
maxIdx  = zeros(nIter,1);

for i = 1:nIter
    [maxComp(i), maxIdx(i)] = max(abs(FList(:,i)));

    fprintf("Iter %2d → max|F| = %.4e at DOF %d\n", ...
        i, maxComp(i), maxIdx(i));
end

%% 4. Detect sudden load redistribution
fprintf("\n---- Load redistribution events ----\n");

for i = 2:nIter
    relChange = norm(FList(:,i) - FList(:,i-1)) / (norm(FList(:,i-1)) + 1e-14);

    if relChange > 0.1
        fprintf("⚠ Iter %2d → large RHS change (%.2f%%)\n", ...
            i, 100*relChange);
    end
end

%% 5. Plot evolution
figure;
plot(Fnorm,'-o');
title('RHS norm evolution');
xlabel('Iteration');
ylabel('||F||');

figure;
plot(maxComp,'-o');
title('Max RHS component evolution');
xlabel('Iteration');
ylabel('max |F_i|');

fprintf("\n==== RHS ANALYSIS END ====\n");