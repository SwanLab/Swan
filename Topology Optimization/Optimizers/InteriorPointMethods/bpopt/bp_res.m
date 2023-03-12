% residual with equality constants and inequality slack variables
classdef bp_res < handle
    properties (Access = public)
        c
    end
    properties (Access = private)
        cRes
        bp
        x
        s
        bL
        bU
    end

    methods (Access = public)
        function obj = bp_res(cParams)
            obj.init(cParams);
        end
        function compute(obj)
            obj.computeBaseResidual();
            obj.computeResidual();
        end
    end
    methods (Access = private)
        function init(obj,cParams)
            obj.bp = cParams.bp;
            obj.x = cParams.x;
            obj.s = cParams.s;
            obj.bL = cParams.bL;
            obj.bU = cParams.bU;
        end
        function computeBaseResidual(obj)
            u.bp = obj.bp;
            u.x = obj.x;
            cResidual = bp_res_stub(u);
            cResidual.create();
            obj.cRes = cResidual.cRes;
        end
        function computeResidual(obj)
            j = 0;
            for i = 1:size(obj.cRes,2)
                if (obj.bU(i) == obj.bL(i))
                    % equality constant
                    obj.cRes(i) = obj.cRes(i) - obj.bL(i);
                else
                    % inequality slack
                    j = j + 1;
                    obj.cRes(i) = obj.cRes(i) - obj.s(j);
                end
            end
            obj.c = obj.cRes;
        end
    end
end