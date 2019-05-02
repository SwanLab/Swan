classdef CC < handle & matlab.mixin.Copyable
    
    properties (Access = public)
        value
        gradient
    end
    
    properties (GetAccess = public, SetAccess = private)
        shapeFunctions
        nSF = 0
    end
    
    properties (Access = private)
        designVariable
        homogVarComputer
    end
    
    methods (Access = public, Abstract)
        updateFields(obj)
    end
    
    methods (Access = public)

        function computeCostAndGradient(obj)
            obj.value = 0;
            obj.gradient = zeros(size(obj.designVariable.value));
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
        
        function obj = init(obj,settings,shapeFuncList,designVariable,homogVarComputer,targetParameters)
            obj.designVariable = designVariable;
            obj.homogVarComputer = homogVarComputer;
            obj.createShapeFunctions(shapeFuncList,settings,targetParameters);
        end        
        
    end
    
    methods (Access = private)
        
        function createShapeFunctions(obj,shapeFunctionNames,settings,targetParameters)
            nShapeFunctions = length(shapeFunctionNames);
            for is = 1:nShapeFunctions
                name          = shapeFunctionNames{is};
                shapeFunction = ShapeFunctional.create(name,settings,obj.designVariable,obj.homogVarComputer,targetParameters);
                obj.append(shapeFunction);
            end
        end
        
        function append(obj,shapeFunction)
            obj.shapeFunctions{obj.nSF+1} = shapeFunction;
            obj.nSF = obj.nSF+1;
        end
              
    end
    
end
