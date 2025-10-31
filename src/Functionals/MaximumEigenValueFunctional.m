classdef MaximumEigenValueFunctional < handle
       
    properties (Access = private)
       density
       eigModes 
       designVariable
       mesh
       filter
       filterAdjoint
       iter
       lambdaOld
       value0
       isCompl
    end
    
    methods (Access = public)
        
        function obj = MaximumEigenValueFunctional(cParams)
            obj.init(cParams)
        end   

        function [f, dfdx] = computeFunctionAndGradient(obj,x) 
            if size(x,2) == 2
                x = x{1};
            end
            xD  = x.obtainDomainFunction();             
            xR = filterDesignVariable(obj,xD{1});
            if obj.isCompl == true
                xR.setFValues(1 - xR.fValues);                             % 1 - F(\chi)
            end
            [lambda,dlambda]= obj.eigModes.computeFunctionAndGradient(xR); % In the 'isCompl' case, the material interpolator recieves 1 - F(\chi)
            f = obj.computeFunction(lambda);
            dfdx{1} = obj.computeGradient(dlambda, lambda);
            if obj.isCompl == true
                dfdx{1}.setFValues(-dfdx{1}.fValues);                      % Chain rule for 1 - F(\chi)
            end
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
            if isfield(cParams,'isCompl')
                obj.isCompl  = cParams.isCompl;
            end
        end

        function J = computeFunction(obj,lambda)
                if isempty(obj.value0)
                    obj.value0  =   1;
                end
                J = -lambda;  
                J = J/obj.value0;
        end

        function dJ = computeGradient(obj, dlambda, lambda)
            if ~isempty(obj.filterAdjoint)
                dJ     = obj.filterAdjoint.compute(dlambda,2);
            else
                dJ        = obj.filter.compute(dlambda,2);
            end
            fValues   = - dJ.fValues;
            dJ.setFValues(fValues/obj.value0);   
        end

        function xR = filterDesignVariable(obj,x)
            xR = obj.filter.compute(x,2);
            if ~isempty(obj.filterAdjoint)
                obj.filterAdjoint.updateFilteredField(obj.filter);
            end
        end

    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = '- First Eigenvalue';
        end
    end  
end