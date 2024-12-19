classdef CGsolver < handle


    properties (Access = private)
      x0  
    end

    methods (Access = public)

        function obj = CGsolver()
            obj.init()
        end

        function x = solve(obj,A,b)
            obj.prepareProblem(b);
            tol = 1e-5;
            maxit = 15000; 
            %L = ichol(A);
            L = diag(sqrt(diag(A)));
            x = pcg(A,b,tol,maxit,L,L',obj.x0);
          %  x = pcg(A,b,tol,maxit,[],[],obj.x0);
            obj.x0 = x;
        end

    end

    methods (Access = private)

        function init(obj)
        end

        function prepareProblem(obj, b)
            n = size(b,1);
            if isempty(obj.x0)
                obj.x0 = zeros(n, 1);
            end        
        end        
    end
end