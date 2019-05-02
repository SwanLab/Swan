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
    end
    
    methods (Access = protected, Abstract)
        updateFields(obj)
    end
    
    methods (Access = public)

        function computeCostAndGradient(obj)
            obj.initValueAndGradient();
            for iSF = 1:length(obj.shapeFunctions)
                obj.shapeFunctions{iSF}.updateTargetParameters();
                obj.shapeFunctions{iSF}.computeCostAndGradient();
                obj.updateFields(iSF);
            end
        end
                       
        function objClone = clone(obj)
            objClone = copy(obj);
        end
        
    end
    
    methods (Access = protected)
        
        function obj = init(obj,cParams)
            obj.nSF   = 0;
            obj.sizeDesigVar = size(cParams.designVariable.value);            
            obj.createShapeFunctions(cParams);
        end        
        
    end
    
    methods (Access = private)
        
        function initValueAndGradient(obj)
            obj.value = 0;
            obj.gradient = zeros(obj.sizeDesigVar);            
        end
        
        function createShapeFunctions(obj,cParams)
            settings         = cParams.settings;
            shapeFuncList    = cParams.shapeFuncList;
            designVariable   = cParams.designVariable;
            homogVarComputer = cParams.homogenizedVarComputer;
            targetParameters = cParams.targetParameters;            
            nShapeFunctions = length(shapeFuncList);
            for is = 1:nShapeFunctions
                name          = shapeFuncList{is};
                shapeFunction = ShapeFunctional.create(name,settings,designVariable,homogVarComputer,targetParameters);
                obj.append(shapeFunction);
            end
        end
        
        function append(obj,shapeFunction)
            obj.shapeFunctions{obj.nSF+1} = shapeFunction;
            obj.nSF = obj.nSF+1;
        end
              
    end
    
end
