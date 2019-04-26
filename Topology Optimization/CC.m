classdef CC < handle
    
    properties (Access = public)
        value
        gradient
        target_parameters
    end
    
    properties (GetAccess = public, SetAccess = private)
        shapeFunctions
        nSF = 0
    end
    
    properties (Access = private)
        designVar
        homogVarComputer
    end
    
    methods (Access = public, Abstract)
        updateFields(obj)
    end
    
    methods (Access = public)

        function computeCostAndGradient(obj,x)
            obj.designVar.value = x;
            obj.value = 0;
            obj.gradient = zeros(size(x));
            for iSF = 1:length(obj.shapeFunctions)
                obj.updateTargetParameters(iSF);
                obj.shapeFunctions{iSF}.computeCostAndGradient();
                obj.updateFields(iSF);
            end
        end
        
    end
    
    methods (Access = protected)
        
        function obj = init(obj,settings,SF_list,designVar,homogVarComputer)
            obj.designVar = designVar;
            obj.homogVarComputer = homogVarComputer;
            obj.createShapeFunctions(SF_list,settings);
        end        
        
    end
    
    methods (Access = private)
        
        function createShapeFunctions(obj,shapeFunctionNames,settings)
            nShapeFunctions = length(shapeFunctionNames);
            for is = 1:nShapeFunctions
                name          = shapeFunctionNames{is};
                shapeFunction = ShapeFunctional.create(name,settings,obj.designVar,obj.homogVarComputer);
                obj.append(shapeFunction);
            end
        end
        
        function append(obj,shapeFunction)
            obj.shapeFunctions{obj.nSF+1} = shapeFunction;
            obj.nSF = obj.nSF+1;
        end
        
        function updateTargetParameters(obj,iSF)
            obj.shapeFunctions{iSF}.target_parameters = obj.target_parameters;
            if contains(class(obj.shapeFunctions{iSF}.filter),'PDE')
                obj.shapeFunctions{iSF}.filter.updateEpsilon(obj.target_parameters.epsilon);
            end
        end
        
    end
    
end
