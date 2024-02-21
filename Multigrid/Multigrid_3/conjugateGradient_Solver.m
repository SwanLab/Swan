function [x,res] = conjugateGradient_Solver(LHS,RHS,x,meshType,maxIter)
    tol = 1e-10;
    %maxIter = 20;
    n = length(RHS);
    %x = RHS; 
    r = RHS - LHS * x; 
    p = r; 
    rsold = r' * r;
    iter = 1;

    hasNotConverged = true;
    
    if strcmp(meshType,'fine')
        
        while iter < maxIter
            Ap = LHS * p;
            alpha = rsold / (p' * Ap);
            x = x + alpha * p;
            r = r - alpha * Ap;
            rsnew = r' * r;

            %hasNotConverged = sqrt(rsnew) > tol;
            hasNotConverged = norm(LHS*x - RHS) > tol;

            p = r + (rsnew / rsold) * p;
            rsold = rsnew;
            iter = iter + 1;
            residu(iter) = norm(LHS*x - RHS); %Ax - b
            res = LHS*x - RHS;

        end
        
    else
        
        while hasNotConverged
            Ap = LHS * p;
            alpha = rsold / (p' * Ap);
            x = x + alpha * p;
            r = r - alpha * Ap;
            rsnew = r' * r;

            %hasNotConverged = sqrt(rsnew) > tol;
            hasNotConverged = norm(LHS*x - RHS) > tol;

            p = r + (rsnew / rsold) * p;
            rsold = rsnew;
            iter = iter + 1;
            residu(iter) = norm(LHS*x - RHS); %Ax - b
            res = LHS*x - RHS;


        end
        
    end
    
end


