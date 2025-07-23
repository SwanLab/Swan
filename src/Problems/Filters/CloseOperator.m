classdef CloseOperator < handle

    properties (Access = private)
        mesh
        dilation
        erotion
    end
    
    methods (Access = public)
        function obj = CloseOperator(cParams)
            obj.init(cParams);
        end
 
        function xFDFE = compute(obj,fun,quadOrder)
            xFD = obj.dilation.compute(fun, quadOrder);
            xFDFE = obj.erotion.compute(xFD, quadOrder);
        end

        function [xF,xFDF] = getFilteredFields(obj)
            xF = obj.dilation.getFilteredField();
            xFDF = obj.erotion.getFilteredField();
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            cParams.filterType = 'FilterAndProject';
            obj.dilation = Filter.create(cParams);
            cParams.eta  = 1 - cParams.eta;
            obj.erotion = Filter.create(cParams);
        end
    end
end