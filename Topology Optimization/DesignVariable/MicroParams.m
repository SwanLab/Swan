classdef MicroParams < DesignVariable
    
    properties (Access = private)
        m1
        m2
    end
    
    methods (Access = public)
        
        function obj = MicroParams(cParams)
            obj.nVariables = 2;
            %obj.m1 = cParams.fun{1};
            %obj.m2 = cParams.fun{2};
            obj.init(cParams);
        end
        
        function update(obj,x)
            x  = obj.splitDesignVariable(x);
            x  = obj.addFixedValues(x);
            x = obj.assambleDesignVariable(x);
            obj.value = x;
        end
        
        function xf = getVariablesToPlot(obj)
            xf    = obj.splitDesignVariable(obj.value);
            xf{3} = obj.computeDensity();
        end
        
        function v = computeVolumeFraction(obj)
            v = obj.computeDensity();
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
            xf = obj.splitDesignVariable(obj.value);
            obj.homogenizedVariablesComputer.computeDensity(xf);
            rho = obj.homogenizedVariablesComputer.rho;
        end        
        
        function createHomogenziedVariableComputer(obj,cParams)
            s.fileName = cParams.homogSettings.fileName;
            s.type = cParams.homogSettings.type;
            s.designVariable = obj;
            h = HomogenizedVarComputer.create(s);
            obj.homogenizedVariablesComputer = h;
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
        
        function x = assambleDesignVariable(obj,xS)
            nVar = obj.nVariables;
            nx = length(xS{1});
            x = zeros(nVar*nx,1);
            for ivar = 1:nVar
                i0 = nx*(ivar-1) + 1;
                iF = nx*ivar;
                xs = xS{ivar};
                x(i0:iF) = xs;
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

