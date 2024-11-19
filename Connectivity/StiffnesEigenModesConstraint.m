classdef StiffnesEigenModesConstraint < handle

    properties (Access = private)
        mesh
        eigenModesFunctional
        shift
        designVariable
        targetEigenValue
    end
    
    methods (Access = public)
        function obj = StiffnesEigenModesConstraint(cParams)
            obj.init(cParams);
            eigen = StiffnessEigenModesComputer(cParams);
            s = cParams;
            s.eigenModes = eigen;
            obj.eigenModesFunctional = MinimumEigenValueFunctional(s);
        end
        
        function [J,dJ] = computeFunctionAndGradient(obj,x)
            [lambda,dlambda] = obj.eigenModesFunctional.computeFunctionAndGradient(x);
            J      = obj.computeFunction(lambda);
            dJ     = obj.computeGradient(dlambda);
        end  

        function updateTargetEigenvalue(obj, new)
            obj.targetEigenValue = new;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh              = cParams.mesh;
            obj.targetEigenValue = cParams.targetEigenValue;
            obj.designVariable    = cParams.designVariable;
            obj.shift             = cParams.shift;
        end

        function J = computeFunction(obj,lambda)
            J    = (obj.shift + obj.targetEigenValue) - lambda;
        end

        function dJ = computeGradient(obj, dlambda)
            dJ = - dlambda;
%             dJ = dJ.project('P1', obj.mesh);
        end

    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Eigenvalue constraint';
        end
    end
end

