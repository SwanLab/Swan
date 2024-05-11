classdef JacobiSolver < handle

    properties (Access = private)
        maxIter
        tol
    end

    methods (Access = public)

        function obj = JacobiSolver(cParams)
            obj.init(cParams);
        end

        function x = solve(obj,LHS,RHS,x)
            n    = length(LHS);
            % x    = zeros(n,1);
            iter = 0;
            D    = diag(diag(LHS));
            T    = LHS - D;
            w    =2/3;

            hasNotConverged = true;

            while iter < obj.maxIter && hasNotConverged
                xold=x;
                x=w*(D\(RHS-T*xold))+(1-w)*xold;
                % if mod(iter,200) == 0
                %     Jacobian_Solver.plotSolution(x,mesh,bc,iter)
                %     res = LHS*x - RHS;
                %     Jacobian_Solver.plotRes(res,mesh,bc,iter)
                % end
                % normVal = norm(x-xold);
                hasNotConverged = norm(LHS*x - RHS) > obj.tol;
                iter  = iter + 1;
                % residu(iter) = normVal;
            end
            % save('residuJacobiZeros.mat','residu')
        end

        % function plotSolution(x,mesh,bc,iter)
        %     xFull = bc.reducedToFullVector(x);
        %     s.fValues = reshape(xFull,2,[])';
        %     s.mesh = mesh;
        %     s.fValues(:,end+1) = 0;
        %     s.ndimf = 3;
        %     xF = P1Function(s);
        %     %xF.plot();
        %     xF.print(['sol',num2str(iter)],'Paraview')
        %     fclose('all');
        % end
        % 
        % function plotRes(res,mesh,bc,numItr)
        %     xFull = bc.reducedToFullVector(res);
        %     s.fValues = reshape(xFull,2,[])';
        %     s.mesh = mesh;
        %     s.fValues(:,end+1) = 0;
        %     s.ndimf = 3;
        %     xF = P1Function(s);
        %     %xF.plot();
        %     xF.print(['Res',num2str(numItr)],'Paraview')
        %     fclose('all');
        % end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.maxIter = cParams.maxIter;
            obj.tol = cParams.tol;
        end
    end
end