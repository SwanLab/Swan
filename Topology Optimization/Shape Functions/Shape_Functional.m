classdef Shape_Functional < handle
    
    properties
        value
        gradient
        target_parameters=struct;
        filter
        Msmooth
        dvolu
        value0
    end
    
    methods (Access = public)
        
        function obj = Shape_Functional()
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,settings)
            obj.filter = Filter.create(settings);
            obj.createMsmoothAndDvolu(settings.filename)
        end
        
        function normalizeFunctionAndGradient(obj)
            obj.normalizeFunctionValue();
            obj.normalizeGradient();
        end
        
    end
    
    methods (Access = private)
        
        function createMsmoothAndDvolu(obj,fileName)
            diffReacProb = DiffReact_Problem(fileName);
            diffReacProb.preProcess;
            obj.Msmooth = diffReacProb.element.M;
            obj.dvolu = diffReacProb.geometry.dvolu;
        end
        
        function normalizeFunctionValue(obj)
            if isempty(obj.value0)
                obj.value0 = obj.value;
            end
            obj.value = obj.value/abs(obj.value0);
        end
        
        function normalizeGradient(obj)
            obj.gradient = obj.gradient/abs(obj.value0);
        end
        
    end
end

