%% BPOPT Solver: Objective function value
classdef bp_obj < handle
    properties (Access = public)
        objectiveFunc
    end
    properties (Access = private)
    end

    methods (Access = public)
        function obj = bp_obj(cParams)
            obj.init(cParams)
        end

        function compute(obj)
            obj.createCases();
        end
    end

    methods (Access = private)
        function init(obj,cParams)
        
        end
        function createCases(obj)
            switch cParams.prob
            case(0)
                y = apm_details(bp.server,bp.app,x);
                ob = y.obj;
                obj.objectiveFunc = ob;
            case(1)
                ob = x(1)*x(4)*(x(1) + x(2) + x(3))  +  x(3);
                obj.objectiveFunc = ob;
            case(2)
                ob = 0;
                n = size(x,2);
                for i = 1:n,
                    ob = ob + (x(i)-i)^2;
                end
                obj.objectiveFunc = ob;
            case(3)
                ob = x(1)^2 - 2*x(1)*x(2) + 4*x(2)^2;
                obj.objectiveFunc = ob;
            case(4)
                ob = x(2)*(5+x(1));
                obj.objectiveFunc = ob;
            case(5)
                ob = x(1)^2-2*x(1)*x(2)+4*x(2)^2;
                obj.objectiveFunc = ob;
            case(6)
                ob = x(1)^2 + 2 * x(2)^2;
                obj.objectiveFunc = ob;
            case(7)
                ob = (x(1)-5)^2 + (x(2)-5)^2;
                obj.objectiveFunc = ob;
            case(8)
                ob = (x(1)-5)^2;
                obj.objectiveFunc = ob;
            case(9)
                ob = (x(1)-5)^2;
                obj.objectiveFunc = ob;
            case(10)
                ob = x(1)^2-2*x(1)*x(2)+4*x(2)^2;
                obj.objectiveFunc = ob;
            end
        end
    end
end
function [ob] = bp_obj(bp,x)

switch bp.prob
    case(0)
        y = apm_details(bp.server,bp.app,x);
        ob = y.obj;
    case(1)
        ob = x(1)*x(4)*(x(1) + x(2) + x(3))  +  x(3);
    case(2)
        ob = 0;
        n = size(x,2);
        for i = 1:n,
            ob = ob + (x(i)-i)^2;
        end
    case(3)
        ob = x(1)^2 - 2*x(1)*x(2) + 4*x(2)^2;
    case(4)
        ob = x(2)*(5+x(1));
    case(5)
        ob = x(1)^2-2*x(1)*x(2)+4*x(2)^2;
    case(6)
        ob = x(1)^2 + 2 * x(2)^2;
    case(7)
        ob = (x(1)-5)^2 + (x(2)-5)^2;
    case(8)
        ob = (x(1)-5)^2;
    case(9)
        ob = (x(1)-5)^2;
    case(10)
        ob = x(1)^2-2*x(1)*x(2)+4*x(2)^2;
end