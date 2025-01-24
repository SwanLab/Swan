classdef ComplianceFunctional < handle

    properties (Access = private)
        value0
    end

    properties (Access = private)
        mesh
        filter
        filterAdjoint
        compliance
        material
        filteredDesignVariable
        gradient
        iter
    end

    methods (Access = public)
        function obj = ComplianceFunctional(cParams)
            obj.init(cParams);
        end

        function [J,dJ] = computeFunctionAndGradient(obj,x)
            iter = x{2};
            x = x{1};
%             
%             if iter > obj.iter
%                 obj.iter = iter;
%                 beta = obj.filter.getBeta();
%                 if iter >= 400 && mod(iter,20)== 0 && beta <= 10
%                     obj.filter.updateBeta(beta + 1.0);
%                     obj.filterAdjoint.updateBeta(beta + 1.0);
%                 end
%             end         

            xD  = x.obtainDomainFunction();
            xR = obj.filterDesignVariable(xD);

            obj.filteredDesignVariable = xR;
            
            obj.material.setDesignVariable(xR);
            [J,dJ] = obj.computeComplianceFunctionAndGradient();
            obj.gradient = dJ;
        end

        function updateFilterParams(obj, beta)
            obj.filter.updateBeta(beta);
            obj.filterAdjoint.updateBeta(beta);
        end

        function beta = getBeta(obj)
            beta = obj.filter.getBeta();
        end

        function x = getDesignVariable(obj)
            x = obj.filteredDesignVariable;
        end

        function x = getGradient(obj)
            x = obj.gradient;
        end
    end
    
    methods (Access = private)
        function init(obj,cParams)
            obj.mesh       = cParams.mesh;
            obj.filter     = cParams.filter;
            obj.material   = cParams.material;
            obj.compliance = cParams.complainceFromConstitutive;
            if isfield(cParams,'value0')
                obj.value0 = cParams.value0;
            end
            if isfield(cParams,'filterAdjoint')
                obj.filterAdjoint = cParams.filterAdjoint;            
            end
            obj.iter = 0;
        end

        function xR = filterDesignVariable(obj,x)
            xR = obj.filter.compute(x,2);
            if ~isempty(obj.filterAdjoint)
                xFiltered = obj.filter.onlyFilter(x,2);
                obj.filterAdjoint.updateFilteredField(xFiltered);
            end
        end

        function [J,dJ] = computeComplianceFunctionAndGradient(obj)
            C   = obj.material.obtainTensor();
            dC  = obj.material.obtainTensorDerivative();
            [J,dJ] = obj.compliance.computeFunctionAndGradient(C,dC);
            if isempty(obj.filterAdjoint) 
                dJ     = obj.filter.compute(dJ,2);
            else
                dJ     = obj.filterAdjoint.compute(dJ,2);
            end

            if isempty(obj.value0)
                obj.value0 = J;
            end
            J          = obj.computeNonDimensionalValue(J);
            dJ.fValues = obj.computeNonDimensionalValue(dJ.fValues);
        end

        function x = computeNonDimensionalValue(obj,x)
            refX = obj.value0;
            x    = x/refX;
        end
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Compliance';
        end
    end
end