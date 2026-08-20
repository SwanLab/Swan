classdef StiffnessEigenModesConstraint < handle

    properties (Access = private)
        mesh
        eigenModesFunctional
        designVariable
        targetEigenValue
        value0
    end
    
    methods (Access = public)
        function obj = StiffnessEigenModesConstraint(cParams)
            obj.init(cParams);
            eigen = StiffnessEigenModesComputer(cParams);
            s = cParams;
            s.eigenModes = eigen;
            obj.eigenModesFunctional = MinimumEigenValueFunctional(s);
        end
        
        function [J,dJ] = computeFunctionAndGradient(obj,x)   
            [lambda,dlambda] = obj.eigenModesFunctional.computeFunctionAndGradient(x);
            J      = obj.computeFunction(lambda);
            dJ{1}     = obj.computeGradient(lambda,dlambda);
        end  
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh              = cParams.mesh;
            obj.targetEigenValue  = cParams.targetEigenValue;
            obj.designVariable    = cParams.designVariable;
        end

        function J = computeFunction(obj,lambda)
              J    = (obj.targetEigenValue - lambda);
              if isempty(obj.value0)
                obj.value0 = 1; 
              end
              J = J/obj.value0;
        end

        function dlambda = computeGradient(obj, lambda, dlambda)
            fValues = - dlambda.fValues/obj.value0;
            dlambda.setFValues(fValues);
        end

    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Eigenvalue constraint';
        end
    end
end

