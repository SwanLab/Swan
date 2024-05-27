classdef MicroParams < DesignVariable

    properties (Access = private)
        density
        structuredMesh
        patchHandle
    end
    
    methods (Access = public)
        
        function obj = MicroParams(cParams)
            obj.nVariables = 2;
            obj.init(cParams);
            obj.density = cParams.density;
            obj.structuredMesh = cParams.structuredMesh;
        end
        
        function update(obj,x)
            x  = obj.splitDesignVariable(x);
            x  = obj.addFixedValues(x);
            for ivar = 1:obj.nVariables                
               obj.fun{ivar}         = obj.fun{ivar}.copy();
               obj.fun{ivar}.fValues = x{ivar};
            end
        end
        
        function xf = getVariablesToPlot(obj)
            xf    = obj.splitDesignVariable(obj.value);
            xf{3} = obj.computeDensity();
        end
        
        function rho = obtainDomainFunction(obj)
            rho = obj.computeDensity();
        end

        function plot(obj)            
            rho   = obj.computeDensity();

            rhoP1D = rho.project('P1D',obj.mesh);
            rhoP1D.plot();    
            colormap(gray)
        end

        
    end
    
    methods (Access = private)
        
        function createValue(obj,cParams)
            [m1,m2] = obj.createM1M2(cParams);
            obj.value = [m1;m2];
        end
        
        
        function createAlpha(obj,cParams)
            if isfield(cParams,'alpha0')
                if ~isempty(cParams.alpha0)
                    obj.alpha = zeros(obj.mesh.ndim,obj.mesh.nelem);
                    obj.alpha(1,:) = cParams.alpha0(1,:);
                    obj.alpha(2,:) = cParams.alpha0(2,:);
                end
            end
        end

        function rho = computeDensity(obj)
          s.operation = @(xV) obj.evaluateDensity(xV);
          rho = DomainFunction(s);
        end
        
        function rhoV = evaluateDensity(obj,xV)
            [mL,cells] = obj.obtainLocalCoord(xV);
            nGaus = size(xV,2);
            rho = obj.density.sampleFunction(mL,cells);
            rhoV(1,:,:) = reshape(rho,nGaus,[]);
        end      

        function [mL,cells] = obtainLocalCoord(obj,xV)
            mx = obj.fun{1};
            my = obj.fun{2};
            mxG = mx.evaluate(xV);
            myG = my.evaluate(xV);
            mG(:,1) = mxG(:);
            mG(:,2) = myG(:);
            [mL,cells] = obj.structuredMesh.obtainLocalFromGlobalCoord(mG);
        end        

        function xS = splitDesignVariable(obj,x)
            nVar = obj.nVariables;
            nx = length(x)/nVar;
            xS = cell(nVar,1);
            for ivar = 1:nVar
                i0 = nx*(ivar-1) + 1;
                iF = nx*ivar;
                xs = x(i0:iF);
                xS{ivar} = xs;
            end
        end
        
        function xVn = addFixedValues(obj,xV)
          xVn = cell(size(xV));
            if ~isempty(obj.isFixed)
                fixedValues = obj.splitDesignVariable(obj.isFixed.values);
                fixNodes = obj.isFixed.nodes;
            else
               fixedValues = cell(numel(xVn));
               fixNodes = [];
            end
            for iVar = 1:obj.nVariables
               xI = xV{iVar};
               fI = fixedValues{iVar};
               xI(fixNodes) = fI(fixNodes);
               xVn{iVar} = xI;
            end
        end
    end
    
end

