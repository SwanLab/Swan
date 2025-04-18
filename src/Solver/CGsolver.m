classdef CGsolver < handle


    properties (Access = private)
      x0  
    end

    methods (Access = public)

        function obj = CGsolver()
            obj.init()
        end

        function x = solve(obj,A,b)
            obj.prepareProblem(A);
            tol = 1e-4;      
            maxit = 15000;

            % Density
            %L = ichol(A);
            %x = pcg(A,b,tol,maxit,L,L',obj.x0);

            % LevelSet de moment..
            x = pcg(A,b,tol,maxit,[],[],obj.x0);
            
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
end