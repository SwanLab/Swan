classdef Cholesky_Direct_Solver < Solver
    methods
        % Analytical Solver (Ax = b)
        function x = solve(obj,A,b)
            A_chol = ichol(A);
            x = obj.preccg(A,b,A_chol',A_chol,10000,1e-6);
        end
    end
    
    methods (Static)

        function x =  preccg(A,b,R1,R2,mmax,tol)
            n = length(A); x = zeros(n,1); r = b-A*x;
            res = zeros(1,mmax); res(1) = dot(r,r); aux = norm(b);
            z = R2\(R1\r);
            z0 = z'*r;
            d = z;
            for m = 1:mmax
                p = A*d;
                xi = z0/dot(d,p);
                x = x+xi*d;
                r = r-xi*p;
                res(m+1) = dot(r,r);
                if (sqrt(res(m+1))<tol*aux)
                    return
                end
                z = R2\(R1\r);
                z1 = z'*r;
                tau = z1/z0;
                d = z+tau*d;
                z0 = z1;
            end
            
            if (m == mmax)
                warning('Maximum number of iterations exceeded.')
            end
        end

    end

end