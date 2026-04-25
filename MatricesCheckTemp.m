clear
close all

%% Load data
S1 = load("LHSListTuron.mat");
S2 = load("LHSListNoTuron.mat");

K_T = S1.LHSList;
K_N = S2.LHSList;

Kanalisis = K_T;

nIter = size(Kanalisis,3);
ndof  = size(Kanalisis,1);

fprintf("==== ANALYSIS START ====\n\n");

%% 1. Rank & condition tracking
fprintf("---- Rank & Condition ----\n");
firstSing = -1;

for i = 1:nIter
    K = Kanalisis(:,:,i);

    r = rank(K);
    c = cond(K);

    fprintf("Iter %2d → rank = %2d, cond = %.2e\n", i, r, c);

    if r < ndof && firstSing < 0
        firstSing = i;
    end
end

if firstSing < 0
    fprintf("\nNo singularity detected.\n");
else
    fprintf("\nFirst singular iteration: %d\n", firstSing);
end

%% 2. Compare with working version
fprintf("\n---- Difference vs NoTuron ----\n");

nCompare = min(size(K_N,3), nIter);

for i = 1:nCompare
    diff = norm(Kanalisis(:,:,i) - K_N(:,:,i), 'fro');
    fprintf("Iter %2d → diff = %.2e\n", i, diff);
end

%% 3. Eigenvalue analysis at failure
if firstSing > 0

    fprintf("\n---- Eigen Analysis (Failing Iteration) ----\n");

    Kbad = Kanalisis(:,:,firstSing);

    [V,D] = eig(Kbad);
    eigvals = diag(D);

    [minEig, idx] = min(abs(eigvals));

    fprintf("Smallest eigenvalue: %.4e\n", minEig);

    mode = V(:,idx);

    %% 4. Plot singular mode
    figure;
    plot(mode,'-o');
    title(['Singular mode at iter ', num2str(firstSing)]);
    xlabel('DOF');
    ylabel('Amplitude');

    %% 5. Compare diagonals
    if firstSing <= size(K_N,3)
        Kgood = K_N(:,:,firstSing);

        figure;
        plot(diag(Kbad),'r-o'); hold on;
        plot(diag(Kgood),'b-o');
        legend('Turon','NoTuron');
        title('Diagonal comparison');
        xlabel('DOF');
        ylabel('Value');
    end
end

%% 6. Minimum eigenvalue evolution
fprintf("\n---- Eigenvalue Evolution ----\n");

minEig = zeros(nIter,1);

for i = 1:nIter
    e = eig(Kanalisis(:,:,i));
    minEig(i) = min(abs(e));
end

figure;
plot(minEig,'-o');
title('Minimum eigenvalue evolution');
xlabel('Iteration');
ylabel('Min |eig|');

fprintf("\n==== ANALYSIS END ====\n");