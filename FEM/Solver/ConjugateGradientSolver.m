classdef ConjugateGradientSolver < handle
    
    properties (Access = private)
        maxIter
        tol
    end
    
    methods (Access = public)
        
        function obj = ConjugateGradientSolver(cParams)
            obj.init(cParams);
        end
        
        function x = solve(obj,LHS,RHS,x)
            r = RHS - LHS * x; 
            p = r; 
            rsold = r' * r;
            iter = 1;

            hasNotConverged = true;

            while iter < obj.maxIter && hasNotConverged
                Ap = LHS * p;
                alpha = rsold / (p' * Ap);
                x = x + alpha * p;
                r = r - alpha * Ap;
                rsnew = r' * r;

                hasNotConverged = norm(LHS*x - RHS) > obj.tol;

                p = r + (rsnew / rsold) * p;
                rsold = rsnew;
                iter = iter + 1;
%                 residu(iter) = norm(LHS*x - RHS); %Ax - b
%                 res = LHS*x - RHS;
                
                %conjugateGradient_Solver.plotSolution(x,mesh,bc,iter)
                
                %conjugateGradient_Solver.plotRes(res,mesh,bc,iter)
            end
            %save('residuConjugateZeros.mat', 'residu')
        end
        
%         function plotSolution(x,mesh,bc,numItr)
%             xFull = bc.reducedToFullVector(x);
%             s.fValues = reshape(xFull,2,[])';
%             s.mesh = mesh;
%             s.fValues(:,end+1) = 0;
%             s.ndimf = 3;
%             xF = P1Function(s);
%             %xF.plot();
%             xF.print(['sol',num2str(numItr)],'Paraview')
%             fclose('all');
%         end
%         
%         function plotRes(res,mesh,bc,numItr)
%             xFull = bc.reducedToFullVector(res);
%             s.fValues = reshape(xFull,2,[])';
%             s.mesh = mesh;
%             s.fValues(:,end+1) = 0;
%             s.ndimf = 3;
%             xF = P1Function(s);
%             %xF.plot();
%             xF.print(['Res',num2str(numItr)],'Paraview')
%             fclose('all');
%         end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.maxIter = cParams.maxIter;
            obj.tol     = cParams.tol;
        end
    end
end