classdef Constraint < handle
    
    properties (Access = public)
        value
        gradient
    end

    properties (Access = private)
        shapeFunctions
     %   Msmooth
    end

    methods (Access = public)
        function obj = Constraint(cParams)
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
            jV  = zeros(nS,1);
            nG  = obj.computeGradientLength(dJc{1});            
            djV = zeros(nG,nS);
            for iS = 1:nS
                jV(iS) = Jc{iS};
                dJs    = dJc{iS};
                dJv = [];                
                for iF = 1:numel(x.fun)
                    dJv = [dJv;dJs{iF}.fValues];                    
                end                
                djV(:,iS) = dJv;
            end
            obj.value    = jV;
%            obj.gradient = obj.Msmooth*djV;
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
                titles{iF} = obj.shapeFunctions{iF}.getTitleToPlot();
            end
        end
    end

    methods (Access = private)
        function obj = init(obj,cParams)
            obj.shapeFunctions = cParams.shapeFunctions;
%            obj.Msmooth        = cParams.Msmooth;
        end
    end
end
