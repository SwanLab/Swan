classdef GMRES_Iterative_Solver < Solver

    properties (Access = private)
        x_Prev
    end
    

    
    methods (Access = public)

        function obj = GMRES_Iterative_Solver(x_Prev)
            obj.init(x_Prev);
        end
        

        %GMRES+Arnoldi Iterative Method (Ax = b) (LHSx=RHS)
        function x_k = solve(obj,LHS,RHS)

            tic
            maxiter = 5000;
            V       = zeros(size(RHS,1),maxiter+1);

            if isempty(obj.x_Prev)
                x_k = RHS;                
            else
                x_k = obj.x_Prev;                
            end

            V(:,1) = x_k/norm(x_k);

            H_n     = zeros(maxiter+1,maxiter);
            z           = LHS*V(:,1);
            H_n(1,1)    = z'*V(:,1);
            v           = z-H_n(1,1)*V(:,1);
            H_n(2,1)    = norm(v);
            V(:,2)      = v/H_n(2,1);
            for j = 2:maxiter
                sum = 0;
                z   = LHS*V(:,j);
                for i = j-1:j
                    H_n(i,j)    = z'*V(:,i);
                    sum         = sum + H_n(i,j)*V(:,i);
                end
                v = z-sum;
                H_n(j+1,j) = norm(v);
                V(:,j+1) = v/H_n(j+1,j);
            end
                e1      = zeros(j+1,1);
                e1(1,1) = 1;      
                C       = H_n;
                D       = norm(RHS)*e1;
                y_k     = C\D;
                x_k     = V(:,1:j)*y_k;
%                 obj.x_k = x_k;                
            toc
        end
    end

    methods (Access = private)

        function init(obj, x_Prev)
            obj.x_Prev = x_Prev;
        end

    end
end