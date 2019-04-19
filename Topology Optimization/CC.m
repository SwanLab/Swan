classdef CC < handle
    
    properties (Access = public)
        value
        gradient
        target_parameters
    end
    
    properties (GetAccess = public, SetAccess = private)
        ShapeFuncs
        nSF = 0
    end
    
    properties (Access = private)
        designVar
    end
    
    methods (Access = public, Abstract)
        updateFields(obj)
    end
    
    methods (Access = public)
        
        function preProcess(obj)
            for iSF = 1:obj.nSF
                obj.ShapeFuncs{iSF}.filter.preProcess;
            end
        end
        
        function computeCostAndGradient(obj, x)
            obj.value = 0;
            obj.gradient = zeros(length(x),1);
            for iSF = 1:length(obj.ShapeFuncs)
                obj.updateTargetParameters(iSF);
                obj.ShapeFuncs{iSF}.computeCostAndGradient(x);
                obj.updateFields(iSF);
            end
        end
        
    end
    
    methods (Access = protected)
        
        function obj = init(obj,settings,SF_list,designVar)
            obj.designVar = designVar;
            obj.createShapeFunctions(SF_list,settings);
        end        
        
    end
    
    methods (Access = private)
        
        function createShapeFunctions(obj,SF_list,settings)
            shapeFactory = ShapeFunctional_Factory();
            for iSF = 1:length(SF_list)
                newShapeFunction = shapeFactory.create(SF_list{iSF},settings,obj.designVar);
                obj.append(newShapeFunction);
            end
        end
        
        function append(obj,shapeFunction)
            obj.ShapeFuncs{obj.nSF+1} = shapeFunction;
            obj.nSF = obj.nSF+1;
        end
        
        function updateTargetParameters(obj,iSF)
            obj.ShapeFuncs{iSF}.target_parameters = obj.target_parameters;
            if contains(class(obj.ShapeFuncs{iSF}.filter),'PDE')
                obj.ShapeFuncs{iSF}.filter.updateEpsilon(obj.target_parameters.epsilon);
            end
        end
        
    end
    
end
