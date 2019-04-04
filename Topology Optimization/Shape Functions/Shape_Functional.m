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
        
        function init(obj,cParams)
            obj.createFilter(cParams.filterParams);
            obj.filter.setupFromGiDFile(cParams.filename,cParams.domainType); 
            obj.createMsmoothAndDvolu(cParams.filename, cParams.domainType);
        end
        
        function normalizeFunctionAndGradient(obj)
            obj.normalizeFunctionValue();
            obj.normalizeGradient();
        end
        
    end
    
    methods (Access = private)
        
        function createFilter(obj,s)
            filterSettings = SettingsFilterFactory(s.filter,s.optimizer);
            obj.filter = FilterFactory().create(filterSettings);
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