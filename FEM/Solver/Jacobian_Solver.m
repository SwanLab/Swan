classdef Jacobian_Solver < Solver

    methods (Static)

        function x = solve(LHS,RHS,mesh,bc)
            normVal = Inf;
            tol = 1e-6;
            n = length(LHS);
            x = zeros(n,1);
            numItr = 0;
            D = diag(diag(LHS));
            T = LHS - D;
            while normVal>tol
            xold=x;
%                 for i=1:n
%                     sigma=0;
%                     for j=1:n
%                         if j~=i
%                             sigma=sigma+LHS(i,j)*x(j);
%                         end
%                     end
                    
            x=D\(RHS-T*x);
            Jacobian_Solver.plotSolution(x,mesh,bc,numItr)
%                 end
            normVal=norm((xold-x)/x);
            numItr = numItr+1;
            end
        
        end
        
        function plotSolution(x,mesh,bc,numItr)
            xFull = bc.reducedToFullVector(x);
            s.fValues = reshape(xFull,2,[])';
            s.mesh = mesh;
            s.fValues(:,end+1) = 0;
            s.ndimf = 3;
            xF = P1Function(s);
            %xF.plot();
            xF.print(['sol',num2str(numItr)],'Paraview')
        end
    end

end