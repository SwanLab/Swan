classdef MinimumEigenValueFunctional < handle
    
    properties (Access = public)
        gradientF
        gradientUN
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
    end
    
    methods (Access = public)
        
        function obj = MinimumEigenValueFunctional(cParams)
            obj.init(cParams)
        end   

        function [lambdas, phis] = computeEigenModes(obj, x, n)
            obj.computeDensity(x);  
            [lambdas, phis] = obj.eigModes.getEigenModesComputer(obj.density,n);  
        end            
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.eigModes       = cParams.eigenModes;
            obj.designVariable = cParams.designVariable;
            obj.mesh           = cParams.mesh;
            if isfield(cParams,'filter')
                obj.filter = cParams.filter;
            end
            if isfield(cParams,'filterAdjoint')
                obj.filterAdjoint  = cParams.filterAdjoint;
            end
            obj.iter = 0;
        end
      
        function computeDensity(obj,x)
            if isempty(obj.filter)
%                 densDomain  = x.fun;
%                 s.operation = @(xV) obj.computeComplementaryDensity(densDomain,xV);
%                 s.mesh = obj.mesh;
%                 densHole = DomainFunction(s);
%                 obj.density = densHole;
            else
                xD  = x.obtainDomainFunction();             % rho
                xR = obj.filterDesignVariable(xD{1});       % FP rho
                xR.setFValues(1 - xR.fValues);              % 1 - FP rho with intermediate densities
                obj.density = xR;
            end

         end
     
        function xR = filterDesignVariable(obj,x)
            if isa(obj.filter, 'HeavisideProjector')
                fValues = obj.filter.project(x);
                xR = FeFunction.create(x.order,fValues,obj.mesh);
            else
                xR = obj.filter.compute(x,2);
            end
            if ~isempty(obj.filterAdjoint)
                xFiltered = obj.filter.getFilteredField();
                obj.filterAdjoint.updateFilteredField(xFiltered);
            end
        end
    end
end