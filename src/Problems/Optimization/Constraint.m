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
                dJc{iF} = obj.mergeGradient(dJ);
            end
            jV  = zeros(nF,1);
            djV = zeros(length(dJc{1}),nF);
            for iF = 1:nF
                jV(iF)    = Jc{iF};
                djV(:,iF) = dJc{iF};
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

        function j = getDesignVariable(obj,i)
            j = obj.shapeFunctions{i}.getDesignVariable();
        end

        function j = getTargetEigenValue(obj,i)
            j = obj.shapeFunctions{i}.getTargetEigenValue();
        end

        function j = getDirichletEigenMode(obj,i)
            j = obj.shapeFunctions{i}.getDirichletEigenMode();
        end

        function j = getGradient(obj,i)
            j = obj.shapeFunctions{i}.getGradient();
        end

        function j = getGradientUN(obj,i)
            j = obj.shapeFunctions{i}.getGradientUN();
        end

        function j = getBeta(obj,i)
            j = obj.shapeFunctions{i}.getBeta();
        end

        function j = getEigenModes(obj)
            j = obj.shapeFunctions{2}.getEigenModes();
        end

        function j = getLambda1(obj)
            j = obj.shapeFunctions{2}.getLambda1();
        end
    end

    methods (Access = private)
        function obj = init(obj,cParams)
            obj.shapeFunctions = cParams.shapeFunctions;
            obj.Msmooth        = cParams.Msmooth;
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