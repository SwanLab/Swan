classdef FilterAndProject < handle

    properties (Access = private)
        filter
        projector
    end

    properties (Access = private)
        mesh
        xFiltered
    end
    
    methods (Access = public)
        function obj = FilterAndProject(cParams)
            obj.init(cParams);
            obj.createFilter(cParams);
            obj.createProjector(cParams);
        end

        function xF = compute(obj,fun,quadOrder)
            obj.xFiltered  = obj.filter.compute(fun,quadOrder);
            xFVal      = obj.projector.project(obj.xFiltered);
            xF         = LagrangianFunction.create(obj.mesh,fun.ndimf,'P1');
            xF.setFValues(xFVal);
        end

        function xF = getFilteredField(obj)
            xF = obj.xFiltered;
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
            if isfield(cParams, 'eta')
                s.eta  = cParams.eta;
            end
            s.beta = cParams.beta;
            obj.projector = HeavisideProjector(s);
        end
    end
end