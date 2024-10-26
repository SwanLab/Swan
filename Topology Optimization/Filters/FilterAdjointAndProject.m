classdef FilterAdjointAndProject < handle

    properties (Access = private)
        filter
        projector
        filteredField
    end

    properties (Access = private)
        mesh
    end
    
    methods (Access = public)
        function obj = FilterAdjointAndProject(cParams)
            obj.init(cParams);
            obj.createFilter(cParams);
            obj.createProjector(cParams);
        end

        function updateFilteredField(obj,xF)
            obj.filteredField = xF;
        end

        function xF = compute(obj,fun,quadOrder)
            sensitVals         = obj.projector.derive(obj.filteredField);
            projSensit         = LagrangianFunction.create(obj.mesh,fun.ndimf,obj.filteredField.order);
            projSensit.fValues = sensitVals;
            regFun             = fun.*projSensit;
            xF                 = obj.filter.compute(regFun,quadOrder);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
        end

        function createFilter(obj,cParams)
            s            = cParams;
            s.filterType = cParams.filterStep;
            obj.filter   = Filter.create(s);
        end

        function createProjector(obj,cParams)
            s.eta  = cParams.eta;
            s.beta = cParams.beta;
            obj.projector = HeavisideProjector(s);
        end
    end
end