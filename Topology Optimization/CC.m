classdef CC < handle 
    
    properties (Access = public)
        value
        gradient
    end
    
    properties (Access = private)
        value0
    end

    properties (Access = private)
        shapeFunctions
        nSF
    end    
    
    methods (Access = protected, Abstract)
        updateFields(obj)
    end
    
    methods (Access = public)

        function computeFunctionAndGradient(obj,x)
            nF = length(obj.shapeFunctions);
            Jc  = cell(nF,1);
            dJc = cell(nF,1);
            for iF = 1:nF
                obj.shapeFunctions{iF}.computeFunctionAndGradient(x);
                Jc{iF}  = obj.weights(iF)*obj.shapeFunctions{iF}.value;
                dJc{iF} = obj.weights(iF)*obj.shapeFunctions{iF}.gradient.fValues;   
            end
            jV  = 0;
            djV = zeros(dJc{1});
            for iF = 1:nF
                jV  = jV  + Jc{iF};
                djV = djV + dJc{iF};
            end
        end
  
    end
    
    methods (Access = protected)
        
        function obj = init(obj,cParams)
            obj.shapeFunctions = cParams.shapeFunctions;
            obj.weights        = cParams.weights;            
            obj.nSF            = length(cParams.shapeFunctions);
        end
        
    end

end
