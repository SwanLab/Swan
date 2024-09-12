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
            nF  = length(obj.shapeFunctions);
            Jc  = cell(nF,1);
            dJc = cell(nF,1);
            for iF = 1:nF
                shI     = obj.shapeFunctions{iF};
                [j,dJ]  = shI.computeFunctionAndGradient(x);
                Jc{iF}  = j;
                dJc{iF} = dJ;
            end
            jV  = zeros(nF,1);
            djV = zeros(length(dJc{1}),nF);
            nF = size(Jc{1,1},2);
            for iF = 1:nF
                jV(iF,1)    = Jc{1,1}(iF);
                djV(:,iF) = dJc{1,1}(:,iF);
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