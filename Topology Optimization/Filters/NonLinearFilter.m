classdef NonLinearFilter < handle
    
    properties (Access = private)
        mesh
        trial
        epsilon
    end

    properties (Access = private)
        % Define new properties here
    end

    methods (Access = public)
        function obj = NonLinearFilter(cParams)
            obj.init(cParams);
            % Construct non-linear stuff...
        end

        function xF = compute(obj,fun,quadOrder)
            xF = LagrangianFunction.create(obj.mesh, 1, obj.trial.order);
            % solve non-linear filter...
            % let's start by creating a factory and define circle case (validation)

            %xF.fValues  = ;
        end

        function obj = updateEpsilon(obj,epsilon)
            if obj.hasEpsilonChanged(epsilon)
                obj.epsilon = epsilon;
                obj.computeLHS();
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.trial   = LagrangianFunction.create(cParams.mesh, 1, cParams.trial.order);
            obj.mesh    = cParams.mesh;
            obj.epsilon = cParams.mesh.computeMeanCellSize();
        end

        function itHas = hasEpsilonChanged(obj,eps)
            var   = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end
    end
end