classdef MinimumEigenValueFunctional < handle
    
    properties (Access = public)
        value
    end
       
    properties (Access = private)
       density
       eigModes 
       designVariable
       mesh
       filter
       filterAdjoint
       iter
       isCompl
    end
    
    methods (Access = public)
        
        function obj = MinimumEigenValueFunctional(cParams)
            obj.init(cParams)
        end   

        function [f, dfdx] = computeFunctionAndGradient(obj,x) 
            %iter = x{2};
            %x = x{1};

            xD  = x.obtainDomainFunction();                  % rho
            xR = obj.filterDesignVariable(xD{1});            % FP rho
            
% 
%             beta = obj.filter.getBeta();
%             if iter >= 100 && mod(iter,100)== 0 && beta <= 10
%                 obj.filter.updateBeta(beta*2.0);
%                 obj.filterAdjoint.updateBeta(beta*2.0);
%             end
% % % 
%             sD.fun      = xR;
%             sD.mesh     = obj.mesh;
%             sD.type     = 'Density';
%             sD.plotting = true;
%             dens        = DesignVariable.create(sD);
%             obj.designVariable = dens;
%             obj.designVariable.plot()

            if obj.isCompl == true
                xR.setFValues(max(min(1 - xR.fValues,1),0)); % 1 - FP
            end
            [f,dfdx]= obj.eigModes.computeFunctionAndGradient(xR);    
            if ~isempty(obj.filterAdjoint)
                dfdx     = obj.filterAdjoint.compute(dfdx,2);
            else
                dfdx     = obj.filter.compute(dfdx,2);
            end
            if obj.isCompl == true
                dfdx.setFValues(-dfdx.fValues);              % Chain rule for (1 - FP)
            end
        end
                
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.eigModes       = cParams.eigenModes;
            obj.designVariable = cParams.designVariable;
            obj.mesh           = cParams.mesh;
            obj.filter = cParams.filter;
            if isfield(cParams,'filterAdjoint')
                obj.filterAdjoint  = cParams.filterAdjoint;
            end
            if isfield(cParams,'isCompl')
                obj.isCompl  = cParams.isCompl;
            end
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
            title = 'MinimumEigenvalue';
        end
    end  
end