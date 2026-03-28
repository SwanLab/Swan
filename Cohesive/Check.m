clear
clc
close all

u     = -0.1;
Kelas = 100;
Kcoh  = 1e8;
Dc    = 0.001;
Df    = 0.1;

R = computeReaction(u, Kelas, Kcoh, Dc, Df)









function R = computeReaction(u, Kelas, Kcoh, Dc, Df)

    fun = @(u2) Kelas*u2 - traction(u - u2, Kcoh, Dc, Df);
    u2 = fzero(fun, u/2);

    Delta = u - u2;
    R = traction(Delta, Kcoh, Dc, Df);

end

function t = traction(Delta, K, Dc, Df)

    if Delta < Dc
        d = 0;
    elseif Delta <= Df
        d = (Delta - Dc)/(Df - Dc);
    else
        d = 1;
    end

    t = (1 - d)*K*Delta;

end