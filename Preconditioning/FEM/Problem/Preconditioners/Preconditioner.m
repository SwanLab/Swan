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
                    M = PreconditionerDirichletNeumann(cParams);
                case {'BlockDiagonal'}
                    M = PreconditionerBlockDiagonal(cParams);
                otherwise
                    error('Invalid preconditioner type.')
            end
        end


        function z = multiplePrec(r,P1,P2,P3,A,b,mesh,bcApplier,uk)


            z1 = P1(r);
            r  = r-A(z1);
            %z2 = P2(r,uk);
            z2 = P2(r);
            r  = r-A(z2);
            z3 = P3(r);
            z  = z1+z2+z3;

            %   J1 = EIFEMtesting.computeTotalEnergy(z1,A,b)
            %   J2 = EIFEMtesting.computeTotalEnergy(z1+z2,A,b)
            %   J3 = EIFEMtesting.computeTotalEnergy(z1+z2+z3,A,b)

            %     EIFEMtesting.plotSolution(z1,mesh,10,10,1,bcApplier,0)
            %     EIFEMtesting.plotSolution(z1+z2,mesh,10,10,2,bcApplier,0)
            %     EIFEMtesting.plotSolution(z1+z2+z3,mesh,10,10,3,bcApplier,0)
            % EIFEMtesting.plotSolution(z,mesh,10,10,3,bcApplier,0)
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


