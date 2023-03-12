% residual with equality constants and inequality slack variables
classdef bp_res < handle
    properties (Access = public)
        c
    end
    properties (Access = private)
        cRes
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
            
        end
        function computeBaseResidual(obj)
            cResidual = bp_res_stub(obj);
            cResidual.create();
            obj.cRes = cResidual.cRes;
        end
        function computeResidual(obj)
            j = 0;
            for i = 1:size(cRes,2),
                if (bU(i)==bL(i)),
                    % equality constant
                    obj.cRes(i) = obj.cRes(i) - bL(i);
                else
                    % inequality slack
                    j = j + 1;
                    obj.cRes(i) = obj.cRes(i) - s(j);
                end
            end
            obj.c = obj.cRes;
        end
    end
end