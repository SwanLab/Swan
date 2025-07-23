classdef StiffnessEigenModesConstraint < handle

    properties (Access = private)
        mesh
        eigenModesFunctional
        designVariable
        targetEigenValue
        filteredDesignVariable
        iter
        filter
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
            iter = x{2};
   
%             if iter > 0 && iter > obj.iter && mod(iter,20)== 0 && obj.targetEigenValue < 2.0
%                 obj.iter = iter;
% %                 obj.targetEigenValue = 2.5;
%                 obj.targetEigenValue = obj.targetEigenValue + 0.1;
%                 disp(obj.targetEigenValue);
%             end
   
            [lambda,dlambda] = obj.eigenModesFunctional.computeFunctionAndGradient(x);
            J      = obj.computeFunction(lambda);
            dJ{1}     = obj.computeGradient(dlambda);
        end  

        function t = getTargetEigenValue(obj)
            t = obj.targetEigenValue;
        end

        function t = getBeta(obj)
            t = obj.eigenModesFunctional.getBeta();
        end

        function updateTargetEigenValue(obj, new)
            obj.targetEigenValue = new;
        end

        function dV = getDesignVariable(obj)
            dV = obj.eigenModesFunctional.getDesignVariable();
        end
       
        function dV = getGradient(obj)
            dV = obj.eigenModesFunctional.getGradient();
        end

        function dV = getGradientUN(obj)
            dV = obj.eigenModesFunctional.getGradientUN();
        end

        function eigenF = getDirichletEigenMode(obj)
            eigenF = obj.eigenModesFunctional.getDirichletEigenMode();
        end

        function eigsCell = getEigenModes(obj)
            eigsCell = obj.eigenModesFunctional.getEigenModes();
        end

        function lambda1 = getLambda1(obj)
            lambda1 = obj.eigenModesFunctional.value;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh              = cParams.mesh;
            obj.targetEigenValue = cParams.targetEigenValue;
            obj.designVariable    = cParams.designVariable;
            obj.iter              = 0.0;
        end

        function J = computeFunction(obj,lambda)
              % if its close to zero, use 1
              J    = (obj.targetEigenValue - lambda);
              if isempty(obj.value0)
                obj.value0 =  1; %abs(J); 
%                 disp(obj.value0)
              end
              J = J/obj.value0;
        end

        function dlambda = computeGradient(obj, dlambda)
            fValues = - dlambda.fValues/obj.value0;
%             dJ      = FeFunction.create(dlambda.order,fValues,obj.mesh);     
            dlambda.setFValues(fValues);
        end

    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Eigenvalue constraint';
        end
    end
end

