classdef Cost < handle 
    
    properties (Access = public)
        value
        gradient
    end

    properties (Access = private)
        shapeFunctions
        weights
        Msmooth
    end

    properties (Access = private)
        shapeValues
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
                [j,dJ]  = shI.computeFunctionAndGradient(x);
                Jc{iF}  = j;
                dJc{iF} = dJ;   
            end
            obj.shapeValues = Jc;
            jV  = 0;
            djV = zeros(size(dJc{1}));
            for iF = 1:nF
                wI  = obj.weights(iF);
                jV  = jV  + wI*Jc{iF};
                djV = djV + wI*dJc{iF};
            end
            obj.value    = jV;
            obj.gradient = obj.Msmooth*djV;
%             obj.gradient = djV;
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
            obj.Msmooth        = cParams.Msmooth;
        end
    end
end