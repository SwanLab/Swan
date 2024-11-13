classdef StiffnesEigenModesConstraint < handle

    properties (Access = private)
        mesh
        minimumEigenValue
        eigenModesFunctional
        shift
        filter
        filterAdjoint
        designVariable
        iter
    end
    
    methods (Access = public)
        function obj = StiffnesEigenModesConstraint(cParams)
            obj.init(cParams);
        end
        
        function [J,dJ] = computeFunctionAndGradient(obj,x)
%             iter = x{2};
%             if iter > 0 && iter > obj.iter && mod(iter,50)== 0 && obj.minimumEigenValue < 0.15
%                 obj.iter = iter;
%                 disp(obj.minimumEigenValue)
%                 obj.minimumEigenValue = obj.minimumEigenValue + 0.01;
%             end
            [lambda,dlambda] = obj.eigenModesFunctional.computeFunctionAndGradient(x);
            J      = obj.computeFunction(lambda);
            dJ     = obj.computeGradient(dlambda);
        end  
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh              = cParams.mesh;
            obj.minimumEigenValue = cParams.minimumEigenValue;
            obj.designVariable    = cParams.designVariable;
            obj.shift             = cParams.shift;
            obj.iter = 0;
            obj.filter            = cParams.filter;
%             s.beta = 20.0;
%             s.eta = 0.5;
%             s.filterStep = 'LUMP';
%             s.mesh = obj.mesh;
%             s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
%             obj.filter = FilterAndProject(s);
%             obj.filterAdjoint = FilterAdjointAndProject(s);

            eigen = StiffnessEigenModesComputer(cParams);
            s.eigenModes = eigen;
            s.designVariable = obj.designVariable;
            s.mesh = obj.mesh;
            obj.eigenModesFunctional = MinimumEigenValueFunctional(s);
        end

        function J = computeFunction(obj,lambda)
            J    = (obj.shift + obj.minimumEigenValue) - lambda;
        end

        function dJ = computeGradient(obj, dlambda)
%             dlambda = obj.filter.compute(dlambda, 2);
%             obj.filterAdjoint.updateFilteredField(dlambda);
%             dJ = obj.filterAdjoint.compute(dlambda,2);

              dJ = obj.filter.compute(dlambda, 2);
        end

    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Eigenvalue constraint';
        end
    end
end

