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
        iter
        filteredDesignVariable
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
%                 if iter >= 400 && mod(iter,20)== 0 && beta <= 40
%                     obj.filter.updateBeta(beta*2.0);
%                     obj.filterAdjoint.updateBeta(beta*2.0);
%                 end
%             end  

            xD  = x.obtainDomainFunction();
            xR = obj.filterFields(xD);
            obj.material.setDesignVariable(xR);
            [J,dJ] = obj.computeComplianceFunctionAndGradient(x);

            obj.filteredDesignVariable = xR{1};
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

        function xR = filterFields(obj,x)
            nDesVar = length(x);
            xR      = cell(nDesVar,1);
            for i = 1:nDesVar
                xR{i} = obj.filter.compute(x{i},2);
                if ~isempty(obj.filterAdjoint)
                    xFiltered = obj.filter.getFilteredField();
                    obj.filterAdjoint.updateFilteredField(xFiltered);
                end
            end
        end

        function [J,dJ] = computeComplianceFunctionAndGradient(obj,x)
            C   = obj.material.obtainTensor();
            dC  = obj.material.obtainTensorDerivative();
            dC  = ChainRule.compute(x,dC);
            [J,dJ] = obj.compliance.computeFunctionAndGradient(C,dC);
            if ~isempty(obj.filterAdjoint)
                dJ{1}     = obj.filterAdjoint.compute(dJ{1},2);
            else
                dJ{1}     = obj.filter.compute(dJ{1},2);
            end
            if isempty(obj.value0)
                obj.value0 = J;
            end
            J  = obj.computeNonDimensionalValue(J);
            dJ = obj.computeNonDimensionalGradient(dJ);
        end

        function x = computeNonDimensionalValue(obj,x)
            refX = obj.value0;
            x    = x/refX;
        end

        function dx = computeNonDimensionalGradient(obj,dx)
            refX = obj.value0;
            for i = 1:length(dx)
                dx{i}.setFValues(dx{i}.fValues/refX);
            end
        end

    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Compliance';
        end
    end
end