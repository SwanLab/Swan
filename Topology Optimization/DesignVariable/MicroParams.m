classdef MicroParams < DesignVariable
    
    properties (Access = private)
        homogenizedVariablesComputer
    end
    
    methods (Access = public)
        
        function obj = MicroParams(cParams)
            obj.nVariables = 2;
            obj.init(cParams);
            obj.createValue(cParams.creatorSettings);
            obj.createAlpha(cParams.creatorSettings)
            obj.createHomogenziedVariableComputer(cParams.creatorSettings);
        end
        
        function update(obj,x)
            %            e
            %
            %              m = obj.mesh;
            % %              %xG = transpose(m.computeBaricenter);
            %              x = m.coord(:,1);
            %              y = m.coord(:,2);
            %
            %             isDirichletPartX = x > -1e-12 & x < 0.1;
            %             isDirichletPartY = y > 0.20 & y < 0.80;
            %             isDirichletPart = isDirichletPartX & isDirichletPartY;
            %             isNeumannPartX = x > (1-0.05) & x < (1+1e-12);
            %             isNeumannPartY = y > 0.40 & y < 0.60;
            %             isNeumannPart = isNeumannPartX & isNeumannPartY;
            %             isForOptimizing = ~isDirichletPart & ~isNeumannPart;
            % %             obj.isElementToOptimize = isForOptimizing;
            xV          = obj.splitDesignVariable(x);
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
            xA = obj.assambleDesignVariable(xVn);
            obj.value = xA;
%             %             xf = obj.getVariablesToPlot();
%             %
%             m1 = xV{1};
%             m2 = xV{2};
%             fixN = [obj.isFixed.nodes obj.isFixed.nodes];
%             m1() = obj.isFixed.values;
%             
%             m1(~isForOptimizing) = 0.011;
%             m2(~isForOptimizing) = 0.011;
            %  obj.value = [m1;m2];
            
            %             s.connec = m.connec(isForOptimizing,:);
            %             s.coord  = m.coord;
            %             m2 = Mesh(s);
            %             m2.plot();
            
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
        
        function [m1,m2] = createM1M2(obj,cParams)
            ndof = length(obj.mesh.coord(:,1));
            m1 = cParams.m1.*ones(ndof,1);
            m2 = cParams.m2.*ones(ndof,1);
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
        
    end
    
end

