classdef GMRES_Iterative_Solver < Solver

    methods (Static)
        %GMRES+Arnoldi Iterative Method (Ax = b) (LHSx=RHS)
        function x_k = solve(LHS,RHS)
            tic
%             TOL     = 1e-3;
            maxiter = 3500;
            V = zeros(size(RHS,1),maxiter+1);
            V(:,1) = RHS/norm(RHS);
            H_n = zeros(maxiter+1,maxiter);

%             res = 100;


            z = LHS*V(:,1);
            H_n(1,1) = z'*V(:,1);
            v = z-H_n(1,1)*V(:,1);
            H_n(2,1) = norm(v);
            V(:,2)   = v/H_n(2,1);

%             while ((res > TOL) && (j <= maxiter))
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

%                 e1 = zeros(j+1,1);
%                 e1(1,1)=1;
%                 H_k = H_n(1:j+1,1:j);
% 
%                 if mod(j,100) == 0
%                 C = H_k'*H_k;
%                 D = H_k'*norm(RHS)*e1;
%                 y_k = C\D;
%                 x_k = V(:,1:j)*y_k;
%                 r_k = RHS-LHS*x_k;
%                 res = norm(r_k)/norm(RHS);
%                 end

%                 j = j+1;
            end
                e1 = zeros(j+1,1);
                e1(1,1)=1;      
                C = H_n;
                D = norm(RHS)*e1;
                y_k = C\D;
                x_k = V(:,1:j)*y_k;

                
%                 load('MATRIU_SVD');                
%                 MATRIU_SVD(:,IT_SVD) = svd(H_n);
%                 IT_SVD = IT_SVD+1; 
%                 save('../GMRES/MATRIU_SVD','MATRIU_SVD','IT_SVD');
            toc
        end
    end
end