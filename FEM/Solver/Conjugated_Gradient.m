classdef Conjugated_Gradient < Solver

    properties (Access = private)
        xOld
        A_trace
    end

    methods (Access = public)
        
        function x = solve(obj,A,b)
            % if isempty(obj.xOld)
            %     obj.A_trace = trace(A);
            %     tol = 0.1; 
            % else
            %     a   = trace(A);
            %     tol = max(1e-8,min(0.1,norm(obj.A_trace - a)));
            %     obj.A_trace = a;
            % end
            % disp('-------TOLERANCE---------')
            % disp(tol)
            % disp('-------------------------')
            if isempty(obj.xOld)
                obj.xOld = zeros(size(A,1),1);
            end
            tol = 1e-6;
            x0 = obj.xOld;
            x  = obj.preccg(A,b,tol,x0);
            obj.xOld = x;
            % obj.incX = x - x0;
        end

    end
    
    methods (Static, Access = private)

        function x =  preccg(A,b,tol,x0)
            % tic
            % L = ichol(A);
            % % L = [];
            % x = pcg(A,b,tol,1e4,L,L',x0);
            % toc
            tic
            x = pcg(A,b,tol,1e4,[],[],x0);
            toc
            tic
            x = A\b;
            toc
        end

    end

end