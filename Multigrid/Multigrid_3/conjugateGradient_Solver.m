function x = conjugateGradient_Solver(LHS,RHS,x)
    tol = 1e-6;
    n = length(RHS);
    %x = RHS; 
    r = RHS - LHS * x; 
    p = r; 
    rsold = r' * r;
    iter = 0;

    hasNotConverged = true;

    while hasNotConverged
        Ap = LHS * p;
        alpha = rsold / (p' * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        rsnew = r' * r;

        %hasNotConverged = sqrt(rsnew) > tol;
        hasNotConverged = max(LHS*x - RHS) > tol;

        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
        iter = iter + 1;
        residu(iter) = norm(LHS*x - RHS); %Ax - b
        res = LHS*x - RHS;
        
        %conjugateGradient_Solver.plotSolution(x,mesh,bc,iter)
        
        %conjugateGradient_Solver.plotRes(res,mesh,bc,iter)
    end
    %save('residuConjugateZeros.mat', 'residu')
end


