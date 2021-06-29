classdef LineSearch < handle
    
    properties (Access = public)
        value
        nTrials
    end
    
    properties (Access = private)
        minValue
        initiator
        designVariable
        objectiveFunction
        initiatorSettings
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = LineSearchFactory();
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = public)
        
        function itIs = isTooSmall(obj)
            itIs = obj.value <= obj.minValue;
        end
        
        function computeStartingValue(obj)
            s = obj.initiatorSettings;
            s.type = 'STANDARD';
            s.designVariable    = obj.designVariable;
            s.objectiveFunction = obj.objectiveFunction;
            s.minValue          = obj.minValue;
            init = LineSearchInitiator.create(s);
            obj.value = init.compute();
            obj.nTrials = 0;            
        end
        
        function computeTrial(obj)
            v = obj.initiator.compute(obj.value);
            obj.value = v;
            obj.nTrials = 0;
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.designVariable    = cParams.designVariable;
            obj.objectiveFunction = cParams.objectiveFunction;
            obj.initiatorSettings = cParams.lineSearchInitiatorSettings;
            obj.minValue = obj.computeMinValue();
            obj.computeLineSearchInitiator();
        end
        
    end
    
    methods (Access = private)
        
        function minV = computeMinValue(obj)
            s = obj.initiatorSettings;
            switch s.optimizerType
                case {'PROJECTED GRADIENT','HAMILTON-JACOBI'}
                    minV = 1e-13;
                case {'SLERP'}
                    minV = 1e-2;
            end                       
        end
        
        function computeLineSearchInitiator(obj)
            s = obj.initiatorSettings;
            s.designVariable    = obj.designVariable;
            s.objectiveFunction = obj.objectiveFunction;
            s.minValue          = obj.minValue;
            obj.initiator = LineSearchInitiator.create(s);
        end
        
        
        
    end
    
    methods (Access = public, Abstract)
        update(obj)
    end
end


