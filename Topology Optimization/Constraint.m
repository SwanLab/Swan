classdef Constraint < handle
    
    properties (Access = public)
        value
        gradient
    end

    properties (Access = private)
        shapeFunctions
        Msmooth
    end

    methods (Access = public)
        function obj = Constraint(cParams)
            obj.init(cParams);
        end

        function computeFunctionAndGradient(obj,x)
            nS  = length(obj.shapeFunctions);
            Jc  = cell(nS,1);
            dJc = cell(nS,1);
            for iF = 1:nS
                shI     = obj.shapeFunctions{iF};
                [j,dJ]  = shI.computeFunctionAndGradient(x);
                Jc{iF}  = j;
                dJc{iF} = dJ;
            end
            jV  = zeros(nS,1);
            djV = zeros(length(dJc{1}),nS);
            for iF = 1:nS
                jV(iF)    = Jc{iF};
                djV(:,iF) = dJc{iF};
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
                titles{iF} = obj.shapeFunctions{iF}.getTitleToPlot();
            end
        end
    end

    methods (Access = private)
        function obj = init(obj,cParams)
            obj.shapeFunctions = cParams.shapeFunctions;
            obj.Msmooth        = cParams.Msmooth;
        end
    end
end