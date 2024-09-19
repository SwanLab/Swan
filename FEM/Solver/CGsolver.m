classdef CGsolver < handle

    % Cuellos de botella:
    % - [CGSolver]: Line 22. I.e Ku=f.
    % - [Projections and PDE Filters] Systems Ax=b with just 1 dof/node may take few time but nothing compared to the elastic problem matrix system.
    % - [ElasticProblem] Elastic stiffness matrix computation
    % - [ElasticProblem] computeStrain and computeStress takes a little bit
    % - [ComplianceFromConstitutive] Line 47 takes few time, i.e, Integrator.compute

    properties (Access = private)
      x0  
    end

    methods (Access = public)

        function obj = CGsolver()
            obj.init()
        end

        function x = solve(obj,A,b)
            obj.prepareProblem(A);
            tol = 1e-5;
            maxit = 15000;
            M = obj.computePreConditionerMatrix(A,'ILU0');
            tic;
            x = pcg(A,b,tol,maxit,M,M',obj.x0); 
            toc;
            obj.x0 = x;
        end

    end

    methods (Access = private)

        function init(obj)
        end

        function prepareProblem(obj, A)
            n = size(A,1);
            if isempty(obj.x0)
                obj.x0 = zeros(n, 1);
            end
        end
    end

    methods (Static, Access = private)
        function M = computePreConditionerMatrix(A,type)
            switch type
                case 'IDENTITY'
                    n = size(A,1);
                    M = sparse(1:n,1:n,ones(1,n),n,n);
                case 'ILU0'
                    M = ichol(A);
            end
        end
    end
end