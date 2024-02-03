function [x,res] = conjugateGradient_Solver(LHS,RHS,x)
    tol = 1e-10;
    maxIter = 20;
    n = length(RHS);
    %x = RHS; 
    r = RHS - LHS * x; 
    p = r; 
    rsold = r' * r;
    iter = 1;

    hasNotConverged = true;

    while iter < maxIter
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
        
        %plotSolution(x,mesh,bc,iter)
        
        %plotRes(res,mesh,bc,iter)
    end
    %save('residuConjugateZeros.mat', 'residu')
end


