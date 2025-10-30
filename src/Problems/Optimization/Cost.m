classdef Cost < handle 
    
    properties (Access = public)
        value
        gradient
    end

    properties (Access = private)
        shapeFunctions
        weights
        Msmooth
        iter
        dofsNonDesign
    end

    properties (Access = private)
        shapeValues
    end

    methods (Access = public)
        function obj = Cost(cParams)
            obj.init(cParams);
            obj.iter = 0.0;
        end

        function computeFunctionAndGradient(obj,x)
%             iter = x{2};
%             if iter > 400 && iter > obj.iter && mod(iter,100) == 0
%                 obj.weights(2) = obj.weights(2)*10;
%                 obj.iter = iter;
%             end

            nF  = length(obj.shapeFunctions);
            Jc  = cell(nF,1);
            dJc = cell(nF,1);
            for iF = 1:nF
                shI     = obj.shapeFunctions{iF};
                [j,dJ]  = shI.computeFunctionAndGradient(x);
                Jc{iF}  = j;
                dJc{iF} = obj.mergeGradient(dJ);   
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

            if ~isempty(obj.dofsNonDesign) 
               obj.gradient(obj.dofsNonDesign) = 0.0;
            end
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
            j = cell2mat(obj.shapeValues(i));
        end

        function j = getDesignVariable(obj,i)
            j = obj.shapeFunctions{i}.getDesignVariable();
        end

        function j = getGradient(obj,i)
            j = obj.shapeFunctions{i}.getGradient();
        end

        function j = getBeta(obj,i)
            j = obj.shapeFunctions{i}.getBeta();
        end
    end
    
    methods (Access = private)
        function init(obj,cParams)
            obj.shapeFunctions = cParams.shapeFunctions;
            obj.weights        = cParams.weights;   
            obj.Msmooth        = cParams.Msmooth;
            if isfield(cParams,'dofsNonDesign') 
                obj.dofsNonDesign = cParams.dofsNonDesign;
            end
        end
    end

    methods (Static,Access=private)
        function dJm = mergeGradient(dJ)
            nDV   = length(dJ);
            nDim1 = length(dJ{1}.fValues);
            dJm   = zeros(nDV*nDim1,1);
            for i = 1:nDV
                ind1           = 1+nDim1*(i-1);
                ind2           = nDim1+nDim1*(i-1);
                indices        = ind1:ind2;
                dJm(indices,1) = dJ{i}.fValues;
            end
        end
    end
end
