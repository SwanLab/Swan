classdef FilterAndProject < handle

    properties (Access = private)
        filter
        projector
    end

    properties (Access = private)
        mesh
    end
    
    methods (Access = public)
        function obj = FilterAndProject(cParams)
            obj.init(cParams);
            obj.createFilter(cParams);
            obj.createProjector(cParams);
        end

        function xF = compute(obj,fun,quadOrder)
            xFiltered  = obj.filter.compute(fun,quadOrder);
            xFVal      = obj.projector.project(xFiltered);
%             xF         = LagrangianFunction.create(obj.mesh,fun.ndimf,fun.order);
            xF         = LagrangianFunction.create(obj.mesh,xFiltered.ndimf,xFiltered.order);
            xF.fValues = xFVal;
        end

        function updateBeta(obj, beta)
            obj.projector.updateBeta(beta);
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