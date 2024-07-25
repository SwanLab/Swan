classdef StiffnesEigenModesConstraint < handle

    properties (Access = private)
        mesh
        minimumEigenValue
        eigenModesFunctional

        designVariable
    end
    
    methods (Access = public)
        function obj = StiffnesEigenModesConstraint(cParams)
            obj.init(cParams);
        end
        
        function [J,dJ] = computeFunctionAndGradient(obj,x)
            [lambda,dlambda] = obj.eigenModesFunctional.computeFunctionAndGradient(x);
            J      = obj.computeFunction(lambda);
            dJ     = obj.computeGradient(dlambda);
        end  
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh             = cParams.mesh;
            obj.minimumEigenValue = cParams.minimumEigenValue;
            obj.designVariable   = cParams.designVariable;

            eigen = StiffnessEigenModesComputer(cParams);
            s.eigenModes = eigen;
            s.designVariable = obj.designVariable;
            obj.eigenModesFunctional = MinimumEigenValueFunctional(s);
             
        end

        function J = computeFunction(obj,lambda)
            J    = obj.minimumEigenValue - lambda;
        end

        function dJ = computeGradient(obj, dlambda)
            dJ      = FeFunction.create('P1', dlambda, obj.mesh);
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Eigenvalue constraint';
        end
    end
end

