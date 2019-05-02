classdef ShapeFunctional < handle
    
    properties
        value
        gradient
        filter
        Msmooth
        dvolu
        value0        
    end
    
    properties (Access = protected)
       homogenizedVariablesComputer 
       nVariables
       designVariable
       target_parameters;                      
    end
        
    methods (Access = public, Static)
        
        function obj = create(type,settings,designVariable,homogVarComputer,targetParameters)
            f = ShapeFunctional_Factory();
            obj = f.create(type,settings,designVariable,homogVarComputer,targetParameters);
        end
        
    end
    
    methods (Access = public)
        
        function updateTargetParameters(obj)
            if contains(class(obj.filter),'PDE')
                obj.filter.updateEpsilon(obj.target_parameters.epsilon);
            end
        end            
                
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.createFilter(cParams.filterParams);
            obj.createMsmoothAndDvolu(cParams.filename, cParams.domainType);
            obj.homogenizedVariablesComputer = cParams.homogVarComputer;
            obj.designVariable = cParams.designVariable;
            obj.target_parameters = cParams.targetParameters;
        end
        
        function normalizeFunctionAndGradient(obj)
            obj.normalizeFunctionValue();
            obj.normalizeGradient();
        end
        
    end
    
    methods (Access = private)
        
        function createFilter(obj,s)
            obj.filter = FilterFactory().create(s);
            obj.filter.preProcess();             
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