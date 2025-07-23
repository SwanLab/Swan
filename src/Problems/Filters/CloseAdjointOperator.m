classdef CloseAdjointOperator < handle

    properties (Access = private)
        dilationAdjoint 
        erotionAdjoint
    end

    properties (Access = private)
        mesh
        xF
        xFDF
        sensDilation
        sensErotion
    end
    
    methods (Access = public)
        function obj = CloseAdjointOperator(cParams)
            obj.init(cParams);
        end
 
        function xEFDF = compute(obj,fun,quadOrder)
            xEF = obj.erotionAdjoint.compute(fun,quadOrder);
            xEFDF = obj.dilationAdjoint.compute(xEF,quadOrder);
        end

        function updateFilteredField(obj,filter)
            [xF,xFDF] = filter.getFilteredFields();
            obj.erotionAdjoint.setFilteredField(xF);
            obj.dilationAdjoint.setFilteredField(xFDF);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            cParams.filterType = 'FilterAdjointAndProject';
            obj.dilationAdjoint = Filter.create(cParams);
            cParams.eta  = 1 - cParams.eta;
            obj.erotionAdjoint = Filter.create(cParams);
        end
    end
end