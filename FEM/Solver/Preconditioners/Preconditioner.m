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
                case {'MODAL'}
                    M = PreconditionerModalApproximation(cParams);
                case {'DirichletNeumann'}
                    M = PreconditionerDirichletNeumann(cParams);
                case {'BNN'}
                    M = 
                otherwise
                    error('Invalid preconditioner type.')
            end
        end


        function z = multiplePrec(r,P1,P2,P3,A)
            z1 = P1(r);
            r  = r-A(z1);
            z2 = P2(r);
            r  = r-A(z2);
            z3 = P3(r);
            z  = z1+z2+z3;
        end

        function z = multiplePrec2(r,P1,P2,A)
            z1 = P1(r);
            r  = r-A(z1);
            z2 = P2(r);
            z  = z1+z2;
        end

        function z = additivePrec(r,P1,P2)
            z1 = P1(r);
            z2 = P2(r);
            z  = z1+z2;
        end

        function x = InexactCG(r,A,P)
            x0 = zeros(size(r));
            factor = 0.99;
            tol = factor*norm(r);
            x = PCG.solve(A,r,x0,P,tol);
            %   tau = 1;

            %tau = @(r,A) 1;
            %   tau = @(r,A) r'*r/(r'*A(r));
            %   x = RichardsonSolver.solve(A,r,x0,P,tol,tau);

        end

    end

end


