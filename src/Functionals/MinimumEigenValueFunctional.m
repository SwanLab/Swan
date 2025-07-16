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

        function [f, dfdx] = computeFunctionAndGradient(obj,x) 
            iter = x{2};
            x = x{1};

%             if iter > obj.iter
%                 obj.iter = iter;
%                 beta = obj.filter.getBeta();
%                 if iter >= 20 && mod(iter,20)== 0 && beta <= 10
%                     obj.filter.updateBeta(beta*2.0);
%                     obj.filterAdjoint.updateBeta(beta*2.0);
%                 end
%             end  

            obj.computeDensity(x);  
            [f,dfdx]= obj.eigModes.computeFunctionAndGradient(obj.density);    
            obj.gradientUN = dfdx;
            if ~isempty(obj.filterAdjoint)
                dfdx     = obj.filterAdjoint.compute(-dfdx,2);
            elseif ~isempty(obj.filter)
                if isa(obj.filter, 'HeavisideProjector')
                    sensitVals         = obj.filter.derive(obj.density);
                    dfdx     = dfdx.project('P1',obj.mesh);
                    dfdx.fValues       = -dfdx.fValues.*sensitVals;
                else
                    dfdx     = obj.filter.compute(-dfdx,2);
%                     dfdx2{1} = dfdx;
                end
            else
                dfdx     = dfdx.project('P1');
            end
            obj.gradientF = dfdx;   
            obj.value = f;
        end

        function [lambdas, phis] = computeEigenModes(obj, x, n)
            obj.computeDensity(x);  
            [lambdas, phis] = obj.eigModes.getEigenModesComputer(obj.density,n);  
        end

        function x = getDesignVariable(obj)
            x = obj.density;
        end

        function x = getBeta(obj)
            x = obj.filter.getBeta();
        end

        function dV = getGradient(obj)
            dV = obj.gradientF;
        end

        function dV = getGradientUN(obj)
            dV = obj.gradientUN;
        end

        function eigenF = getDirichletEigenMode(obj)
            eigenF = obj.eigModes.phiDCont;
        end

        function eigsF = getEigenModes(obj)
            eigsF = obj.eigModes.eigsF;
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
%                 xR.setFValues(1 - max(0,min(1,round(xR.fValues)))); % 1 - FP rho without intermediate densities
                obj.density = xR;
            end
%             s.fun  = obj.density;
%             s.mesh = obj.mesh;
%             s.type = 'Density';
%             s.plotting = true;
%             dens    = DesignVariable.create(s);
%             dens.plot();
         end
      
        function rho = computeComplementaryDensity(obj,fun,xV)
            rho = fun.evaluate(xV);
            rho = 1 - rho;
            rho = round(rho);
            rho = max(0,min(1,rho));
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
    
    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'MinimumEigenvalue';
        end
    end  
end