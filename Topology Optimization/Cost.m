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
            nS  = length(obj.shapeFunctions);
            Jc  = cell(nS,1);
            dJc = cell(nS,1);
            for iS = 1:nS
                shI     = obj.shapeFunctions{iS};
                [j,dJ]  = shI.computeFunctionAndGradient(x);
                Jc{iS}  = j;
                dJc{iS} = dJ;   
            end
            obj.shapeValues = Jc;
            jV  = 0;
            nG  = obj.computeGradientLength(dJc{1});
            djV = zeros(nG,1);            
            for iS = 1:nS
                wI  = obj.weights(iS);
                jV  = jV  + wI*Jc{iS};
                dJs = dJc{iS};
                dJv = [];                
                for iF = 1:numel(x.fun)
                    dJv = [dJv;dJs{iF}.fValues];                    
                end
                djV = djV + wI*dJv;
            end
            obj.value    = jV;
            obj.gradient = djV;
        end

        function nG = computeGradientLength(obj,dJ)
            nF = numel(dJ);
            nG = 0;
            for iF = 1:nF
                nG = nG + length(dJ{iF}.fValues);
            end
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
