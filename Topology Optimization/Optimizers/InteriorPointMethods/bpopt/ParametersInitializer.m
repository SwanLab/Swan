%% BPOPT Solver: Initialize variables and equations
%
%% Problem 0
%     APMonitor Modeling Langauge (http://apmonitor.com)
%% Problem 1
%     min   x1*x4*(x1 + x2 + x3)  +  x3
%     s.t.  x1*x2*x3*x4                   >=  25
%           x1^2 + x2^2 + x3^2 + x4^2  =  40
%           1 <=  x1,x2,x3,x4  <= 5
%
%     Starting point:
%        x = (1, 5, 5, 1)
%
%     Optimal solution:
%        x = (1.00000000, 4.74299963, 3.82114998, 1.37940829)
%% Problem 2
%     min   sum (x(i)-i)^2
%     s.t.  sum(x(i)) = 10
%
%     Starting point:
%        x = zeros
%% Problem 3
%     min   x1^2 - 2*x1*x2 + 4*x2^2
%     s.t.  0.1*x1 - x2 > 1
%% Problem 4
%     min   x2*(5+x1)
%     s.t.  x1 * x2 >= 5
%           x1^2 + x2^2 <=20
%% Problem 5
%     min   x1^2 - 2*x1*x2 + 4*x2^2
%     s.t.  0.1*x1 - x2 > 1
%           -10*x(1)+x(2) > 1
%% Problem 6
%     min   x1^2 + (2 * x2)^2
%     s.t.  2 * x1 + x2 <= 9
%           x1 + 2 * x2  = 10
%% Problem 7
%     min   (x1-5)^2 + (x2-5)^2
%     s.t.  x1 + x2 = 10
%           0 <= x1 <= 10
%           0 <= x2 <= 9
%% Problem 8 
%     min   (x1-5)^2
%     s.t.  9 - x1^2 - x2 = 0
%           0 <= x1 <= 10
%           x2 >= 0
%% Problem 9
%     min   (x1-5)^2
%     s.t.  x1^2 <= 9
%           0 <= x1 <= 10
%% Problem 10
%     min   x1^2 - 2*x1*x2 + 4*x2^2
%     s.t.  0.1 * x1 - x2 > 1
classdef ParametersInitializer < handle
    properties (Access = public)
        x0C
        xLC
        xUC
        bLC
        bUC
    end
    properties (Access = private)
        bp
    end

    methods (Access = public)
        function obj = ParametersInitializer(cParams)
            obj.init(cParams);
        end
        function create(obj)
            obj.selectProblem();
        end
    end
    methods (Access = private)
        function init(obj,cParams)
            obj.bp = cParams;
        end
        function selectProblem(obj)
            switch obj.bp.prob
            case(0)
                %y = apm_details(obj.bp.server,obj.bp.app);
                %x0 = y.var.val';
                %xL = y.var.lb';
                %xU = y.var.ub';
                %bL = y.eqn.lb';
                %bU = y.eqn.ub';
            case(1)
                x0 = [1, 5, 5, 1];
                xL = ones(1,4);
                xU = 5*ones(1,4);
                bL = [25 40];
                bU = [1e10 40];
            case(2)
                n = 2;
                x0 = ones(1,n);
                xL = 0*ones(1,n);
                xU = 100*ones(1,n);
                xU(1) = 10;
                bL = 10;
                bU = 10;
            case(3)
                x0 = [5, -3];
                xL = [-10 -10];
                xU = [10 10];
                bL = 1;
                bU = 1e20;
            case(4)
                x0 = [2, 3];
                xL = [1, 1];
                xU = [5, 5];
                bL = [5, -1e20];
                bU = [1e20, 20];
            case(5)
                x0 = [-7 -5];
                xL = [-1e20 -1e20];
                xU = [1e20 1e20];
                bL = [1 1];
                bU = [1e20 1e20];
            case(6)
                x0 = [3 2];
                xL = [0 0];
                xU = [1e20 1e20];
                bL = [-1e20 10];
                bU = [9 10];
            case(7)
                x0 = [1 1];
                xL = [0 0];
                xU = [10 9];
                bL = 10;
                bU = 10;
            case(8)
                x0 = [4 4];
                xL = [0 0];
                xU = [10 1e10];
                bL = 0;
                bU = 0;
            case(9)
                x0 = 4;
                xL = 0;
                xU = 10;
                bL = -1e10;
                bU = 9;
            case(10)
                x0 = [-5 -5];
                xL = [-1e20 -1e20];
                xU = [1e20 1e20];
                bL = [1];
                bU = [1e20];
            end
            obj.x0C = x0;
            obj.xLC = xL;
            obj.xUC = xU;
            obj.bLC = bL;
            obj.bUC = bU;
        end
    end
end