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
            obj.createFilter(settings);
            obj.createMsmoothAndDvolu(settings.filename,settings.ptype);
        end
        
        function normalizeFunctionAndGradient(obj)
            obj.normalizeFunctionValue();
            obj.normalizeGradient();
        end
        
    end
    
    methods (Access = private)
        
        function createFilter(obj,settings)
            obj.filter = FilterFactory().create(settings.filter,settings.optimizer);
            obj.filter.setupFromGiDFile(settings.filename,settings.ptype); 
        end
        
        function createMsmoothAndDvolu(obj,fileName,scale)
            diffReacProb = obj.createDiffReactProb(scale);
            diffReacProb.setupFromGiDFile(fileName);
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
    
    methods (Static, Access = private)
        function diffReacProb = createDiffReactProb(scale)
            switch scale
                case 'MACRO'
                    diffReacProb = DiffReact_Problem;
                case 'MICRO'
                    diffReacProb = DiffReact_Problem_Micro;
            end
        end
    end
end