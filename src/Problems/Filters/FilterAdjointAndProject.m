classdef FilterAdjointAndProject < handle

    properties (Access = private)
        filter
        projector
        filteredField
        projSensit
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

        function xF = compute(obj,fun,quadOrder)
            if isempty(obj.projSensit)
                obj.projSensit = LagrangianFunction.create(obj.mesh,fun.ndimf,'P1');
            end
            sensitVals = obj.projector.derive(obj.filteredField);
            obj.projSensit.setFValues(sensitVals);
            regFun = fun.*obj.projSensit;
            xF     = obj.filter.compute(regFun,quadOrder);
        end

%         function updateFilteredField(obj,xF)
%             obj.filteredField = xF;
%         end

        function updateFilteredField(obj,filter)
           xF = filter.getFilteredField();
           obj.setFilteredField(xF);
        end

        function setFilteredField(obj,xF)
            obj.filteredField = xF;
        end

        function updateBeta(obj, beta)
            obj.projector.updateBeta(beta);
        end

        function updateEpsilon(obj,epsilon)
            obj.filter.updateEpsilon(epsilon);
        end

        function beta = getBeta(obj)
            beta = obj.projector.getBeta();
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
            s.beta = cParams.beta;
            if isfield(cParams, 'eta')
                s.eta  = cParams.eta;
            end
            obj.projector = HeavisideProjector(s);
        end

    end
end