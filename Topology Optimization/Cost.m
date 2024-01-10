classdef Cost < handle 
    
    properties (Access = public)
        value
        gradient
    end

    properties (Access = private)
        shapeFunctions
        weights
    end    
    
    methods (Access = public)

        function computeFunctionAndGradient(obj,x)
            nF  = length(cParams.shapeFunctions);
            Jc  = cell(nF,1);
            dJc = cell(nF,1);
            for iF = 1:nF
                shI     = obj.shapeFunctions{iF};
                wI      = obj.weights(iF);
                [j,dJ]  = shI.computeFunctionAndGradient(x);
                Jc{iF}  = wI*j;
                dJc{iF} = wI*dJ.fvalues;   
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
        end
        
    end

end
