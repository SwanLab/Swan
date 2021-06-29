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
            obj.sizeDesigVar = size(cParams.designVar.value);            
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
                s.homogVarComputer = cParams.homogenizedVarComputer;
                s.targetParameters = cParams.targetParameters;
                shapeFunction = ShapeFunctional.create(s);
                obj.append(shapeFunction);
            end
        end
        
        
        function append(obj,shapeFunction)
            obj.shapeFunctions{obj.nSF+1} = shapeFunction;
            obj.nSF = obj.nSF+1;
        end
              
    end
    
end
