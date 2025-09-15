classdef Preconditioner < handle

    methods (Access = public, Static)

        function M = create(cParams)

            switch cParams.type
                case {'Jacobi'}
                    M = PreconditionerJacobi(cParams);
                case {'GaussSeidel'}
                    M = PreconditionerGaussSeidel(cParams);
                case {'ILU'}
                    M = PreconditionerILU(cParams);
                case {'EIFEM'}
                    M = PreconditionerEIFEM(cParams);
                case {'EIFEMcont'}
                    M = PreconditionerEIFEMcontinous(cParams);
                case {'MODAL'}
                    M = PreconditionerModalApproximation(cParams);
                case {'DirichletNeumann'}
                    M = PreconditionerDirichletNeumann(cParams);
                case {'BNN'}
                    M = PreconditionerDirichletNeumann(cParams);
                case {'BlockDiagonal'}
                    M = PreconditionerBlockDiagonal(cParams);
                otherwise
                    error('Invalid preconditioner type.')
            end
        end

        function z = multiplePrec(r, A, varargin)
            z = 0;
            for k = 1:length(varargin)
                Pk = varargin{k};
                zk = Pk(r);
                r  = r - A(zk);
                z  = z + zk;
            end
        end


        % function z = multiplePrec(r,P1,P2,P3,A)
        %     z1 = P1(r);
        %     r  = r-A(z1);
        %     z2 = P2(r);
        %     r  = r-A(z2);
        %     z3 = P3(r);
        %     z  = z1+z2+z3;
        % end
        % 
        % function z = multiplePrec2(r,P1,P2,A)
        %     z1 = P1(r);
        %     r  = r-A(z1);
        %     z2 = P2(r);
        %     z  = z1+z2;
        % end

        function z = additivePrec(r,P1,P2)
            z1 = P1(r);
            z2 = P2(r);
            z  = z1+z2;
        end

        function x = InexactCG(r,A,P,b)
            x0 = zeros(size(r));
           
            factor = 1-1e-4;%0.99;
            tol = factor*norm(r);
%             tol = 1;
            
      %      x = PCG.solve(A,r,x0,P,tol);
            
            %tau = @(r,A) 1;
           %tau = @(z,r,A) r'*z/(z'*A(z)); 
           tau = @(z,r,A) 1; 
    %       tau = @(z,r,A) r'*r/(r'*A(r));           
           x = RichardsonSolver.solve(A,r,x0,P,tol,tau);
        end

    end

end


