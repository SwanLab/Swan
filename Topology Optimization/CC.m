classdef CC < handle & matlab.mixin.Copyable
    
    properties (Access = public)
        value
        gradient
    end
    
    properties (GetAccess = public, SetAccess = private)
        shapeFunctions
        nSF
    end
    
    properties (Access = private)
        sizeDesigVar
        valueOld
        gradientOld
    end
    
    methods (Access = protected, Abstract)
        updateFields(obj)
    end
    
    methods (Access = public)

        function computeFunctionAndGradient(obj)
            obj.initValueAndGradient();
            for iSF = 1:length(obj.shapeFunctions)
                obj.shapeFunctions{iSF}.updateTargetParameters();
                obj.shapeFunctions{iSF}.computeFunctionAndGradient();
                obj.updateFields(iSF);
            end
        end
        
        function computeFunction(obj)
            obj.initValueAndGradient();
            for iSF = 1:length(obj.shapeFunctions)
                obj.shapeFunctions{iSF}.updateTargetParameters();
                obj.shapeFunctions{iSF}.computeFunction();
                obj.updateFields(iSF);
            end
        end
        
        function objClone = clone(obj)
            objClone = copy(obj);
        end
        
        function restart(obj)
            obj.value    = obj.valueOld;
            obj.gradient = obj.gradientOld;
        end
        
        function updateOld(obj)
           obj.valueOld    = obj.value;
           obj.gradientOld = obj.gradient;
        end
        
    end
    
    methods (Access = protected)
        
        function obj = init(obj,cParams)
            obj.nSF   = 0;
            obj.sizeDesigVar = size(cParams.designVar.fun.fValues);
            obj.createShapeFunctions(cParams);
        end
        
    end
    
    methods (Access = private)
        
        function initValueAndGradient(obj)
            obj.value = 0;
            obj.gradient = zeros(obj.sizeDesigVar);
        end
        
        function createShapeFunctions(obj,cParams)
            nS = cParams.nShapeFuncs;
            for iS = 1:nS
                s = cParams.shapeFuncSettings{iS};
                s.designVariable = cParams.designVar;
                s.femSettings.physicalProblem = cParams.physicalProblem;
                s.shNumber = iS;

                if isprop(cParams,'homogenizedVarComputer')
                    s.homogVarComputer = cParams.homogenizedVarComputer;
                end
                if isprop(cParams,'targetParameters')
                    s.targetParameters = cParams.targetParameters;
                end
                s.femSettings.designVariableFilter = obj.createFilter(s);
                s.femSettings.gradientFilter       = obj.createFilter(s);
                shapeFunction = ShapeFunctional.create(s);
                obj.append(shapeFunction);
            end
        end

        function append(obj,shapeFunction)
            obj.shapeFunctions{obj.nSF+1} = shapeFunction;
            obj.nSF = obj.nSF+1;
        end
       
    end

    methods (Static, Access = private)
        function filter = createFilter(cParams)
            functionalType = cParams.type;
            s              = cParams.filterParams.femSettings;
            switch functionalType
                case {'perimeter','perimeterConstraint'}
                    s.filterType   = 'PDE';
                    s.boundaryType = 'Robin';
                case {'perimeterInterior'}
                    s.filterType = 'PDE';
                    if s.scale == "MICRO"
                        s.boundaryType = 'Periodic';
                    else
                        s.boundaryType = 'Neumann';
                    end
                case {'anisotropicPerimeter2D'}
                    s.filterType        = 'PDE';
                    s.boundaryType      = 'Robin';
                    u                   = 45;
                    s.CAnisotropic      = [tand(u),0;0,1/tand(u)];
                    s.aniAlphaDeg       = 90;
                    s.metric            = 'Anisotropy';
                case {'anisotropicPerimeterInterior2D'}
                    s.filterType        = 'PDE';
                    u                   = 45;
                    s.CAnisotropic      = [tand(u),0;0,1/tand(u)];
                    s.aniAlphaDeg       = 90;
                    s.metric            = 'Anisotropy';
                    if s.scale == "MICRO"
                        s.boundaryType = 'Periodic';
                    else
                        s.boundaryType = 'Neumann';
                    end
                otherwise
                    s.filterType = cParams.filterParams.filterType;
            end
            s.mesh  = cParams.designVariable.mesh;
            s.test  = P0Function.create(s.mesh,1);
            s.trial = P1Function.create(s.mesh,1);
            filter  = Filter.create(s);
        end
    end

end
