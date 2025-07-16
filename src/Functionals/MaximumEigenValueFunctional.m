classdef MaximumEigenValueFunctional < handle
       
    properties (Access = private)
       density
       eigModes 
       designVariable
       mesh
       filter
       filterAdjoint
       iter
       value0
    end
    
    methods (Access = public)
        
        function obj = MaximumEigenValueFunctional(cParams)
            obj.init(cParams)
        end   

        function [f, dfdx] = computeFunctionAndGradient(obj,x) 
            x = x{1};
            xD  = x.obtainDomainFunction();             
            xR = obj.filter.compute(xD{1},2);   
            if ~isempty(obj.filterAdjoint)
                xFiltered = obj.filter.getFilteredField();
                obj.filterAdjoint.updateFilteredField(xFiltered);
            end
            [lambda,dlambda]= obj.eigModes.computeFunctionAndGradient(xR);
            f = obj.computeFunction(lambda);
            dfdx{1} = obj.computeGradient(dlambda);
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.eigModes       = StiffnessEigenModesComputer(cParams);
            obj.designVariable = cParams.designVariable;
            obj.mesh           = cParams.mesh;
            obj.filter = cParams.filter;
            if isfield(cParams,'filterAdjoint')
                obj.filterAdjoint = cParams.filterAdjoint;
            end
        end

        function J = computeFunction(obj,lambda)
%                  J = 1/lambda;
                J = -lambda;
                if isempty(obj.value0)
                   obj.value0 = abs(J);
                end
                J = J/obj.value0;
        end

        function dJ = computeGradient(obj, dlambda)
            if ~isempty(obj.filterAdjoint)
                dJ     = obj.filterAdjoint.compute(dlambda,2);
            else
                dJ        = obj.filter.compute(dlambda,2);
            end
            fValues   = - dJ.fValues;
%             fValues   =  -1./dJ.fValues.^2; %
            dJ.setFValues(fValues/obj.value0);   
        end

    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = '- First Eigenvalue';
        end
    end  
end