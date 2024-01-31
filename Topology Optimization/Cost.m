classdef Cost < handle 
    
    properties (Access = public)
        value
        gradient
    end

    properties (Access = private)
        shapeFunctions
        shapeValues
        weights
    end

    methods (Access = public)
        function obj = Cost(cParams)
            obj.init(cParams);
        end

        function computeFunctionAndGradient(obj,x)
            nF  = length(obj.shapeFunctions);
            Jc  = cell(nF,1);
            dJc = cell(nF,1);
            for iF = 1:nF
                shI     = obj.shapeFunctions{iF};
                wI      = obj.weights(iF);
                [j,dJ]  = shI.computeFunctionAndGradient(x);
                obj.shapeValues{iF} = j;
                Jc{iF}  = wI*j;
                dJc{iF} = wI*dJ.fValues;   
            end
            jV  = 0;
            djV = zeros(size(dJc{1}));
            for iF = 1:nF
                jV  = jV  + Jc{iF};
                djV = djV + dJc{iF};
            end
            obj.value    = jV;
            obj.gradient = djV;
        end

        function nF = obtainNumberFields(obj)
            nF = length(obj.shapeFunctions);
        end

        function titles = getTitleFields(obj)
            nF = length(obj.shapeFunctions);
            titles = cell(nF,1);
            for iF = 1:nF
                wI         = obj.weights(iF);
                titleF     = obj.shapeFunctions{iF}.getTitleToPlot();
                titles{iF} = [titleF,' (w=',int2str(wI),')'];
            end
        end

        function j = getFields(obj,i)
            j = obj.shapeValues{i};
        end
    end
    
    methods (Access = private)
        function obj = init(obj,cParams)
            obj.shapeFunctions = cParams.shapeFunctions;
            obj.weights        = cParams.weights;            
        end
    end
end
