function [x,res,residuFine,iterRes] = conjugateGradient_Solver(LHS,RHS,x,malla,residuFine,iterRes,maxIter)
    tol = 1e-10;
    %maxIter = 20;
    n = length(RHS);
    %x = RHS; 
    r = RHS - LHS * x; 
    p = r; 
    rsold = r' * r;
    iter = 1;

    hasNotConverged = true;
    
    if strcmp(malla,'fine')
        
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

            if strcmp(malla,'fine')
                residuFine(iterRes) = norm(LHS*x - RHS);
                iterRes = iterRes + 1;
            end

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

            if strcmp(malla,'fine')
                residuFine(iterRes) = norm(LHS*x - RHS);
                iterRes = iterRes + 1;
            end

        end
        
    end
    
end


