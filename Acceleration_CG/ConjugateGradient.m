classdef ConjugateGradient < handle

    properties (Access = private)
        xOld
        tol
        maxIters
        tolFactor
        xBest
        volume
        xConverged
    end

    methods (Access = public)

        function obj = ConjugateGradient(cParams)
            obj.tol      = cParams.tol;
            obj.volume   = cParams.volume;
            obj.maxIters = 6e3; % - Modifiable
        end

        function x = solve(obj,A,b)
                tic
            if isempty(obj.xOld) % - Might be replaced in the future by an stored initial guess
                x         = A\b;
                obj.xBest = x;
                obj.xConverged = x;
                disp('DIRECT METHOD')
            else
                t = obj.tol.val;
                x = obj.xOld;
                rBest = norm(b - A*obj.xBest);
                rOld  = norm(b - A*obj.xOld);
                if rBest <= rOld
                    x = obj.xBest;
                else
                    x = obj.xOld;
                    obj.xBest = obj.xOld;
                end
                if norm(b - A*obj.xBest) > t && norm(b - A*obj.xOld) > t %norm(b-A*x)/norm(b) > t %
                    % x = obj.xConverged;
                    % Fitting params
                    thetay = [8374.61988218320;
                             -14755.8486727114;
                              7173.21553965786];
                    % - 
                    thetax = [1177.73782766053;
                              -2144.56933446969;
                              1210.63290534893];
                    % Obtaining factors
                    v    = obj.volume.currentVal;
                    yMax = [v^2,v,1]*thetay;
                    xMax = [v^2,v,1]*thetax;

                    % Solving the problem
                    % [xTrial,~,~,it_1] = pcg(A,b,t,100,[],[],x);
                    % Scaling
                    alphay = yMax/max(abs(x(2:2:end)));
                    alphax = xMax/max(abs(x(1:2:end)));
                    x(2:2:end) = alphay*x(2:2:end);
                    x(1:2:end) = alphax*x(1:2:end);
                    [x,~,~,it] = pcg(A,b,t,obj.maxIters,[],[],x);
                else
                    it = 0;
                end
                disp('Iter: ' + string(it))
            end
            obj.xOld = x;
            disp('Tol: ' + string(obj.tol.val))
            disp('Convergence time: ' + string(toc) + ' s')
            disp('------------------')
        end

        % function x = solve(obj,A,b)
        %     tic
        %     if isempty(obj.xOld) % - Might be replaced in the future by an stored initial guess
        %         % obj.defineCorrectionFactor(b);
        %         % x         = A\b;
        %         % obj.xBest = x;
        %         % disp('DIRECT METHOD')
        %         x = zeros(numel(b),1);
        %         [x,~,~,it] = pcg(A,b,obj.tol.val,obj.maxIters,[],[],x);
        %         obj.xBest = x;
        %     else
        %         if obj.tol.val > 0
        %             t = obj.tol.val;%*obj.tolFactor;
        %             rBest = norm(b - A*obj.xBest);
        %             rOld  = norm(b - A*obj.xOld);
        %             if rBest <= rOld
        %                 x = obj.xBest;
        %             else
        %                 x = obj.xOld;
        %                 obj.xBest = obj.xOld;
        %             end
        %             [x,~,~,it] = pcg(A,b,t,obj.maxIters,[],[],x);
        %             disp('Iter: ' + string(it))
        %         else
        %             x = A\b;
        %         end
        %     end
        %     obj.xOld = x;
        %     disp('Tol: ' + string(obj.tol.val))
        %     disp('Convergence time: ' + string(toc) + ' s')
        %     disp('------------------')
        % end

    end

    methods (Access = private)

        function defineCorrectionFactor(obj,b)
            obj.tolFactor = 1;
        end

    end
end