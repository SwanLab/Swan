classdef MultimaterialGradientComputer < handle

    properties (Access = private)
        designVariable
    end

    methods (Access = public)
        function obj = MultimaterialGradientComputer(desVar)
            obj.init(desVar);
        end

        function dt = compute(obj,TD)
            x      = obj.designVariable;
            tfi    = x.obtainGlobalDomainFunction();
            tfiDer = x.obtainDomainFunctionDerivatives();
            tfi    = obj.expandVariable(tfi,TD{1,1}); % Provisional until ndimf vector
            tfiDer = obj.expandVariable(tfiDer,TD{1,1}); % Provisional until ndimf vector
            nLS    = length(tfi)-1;
            dt     = cell(nLS,1);
            for k = 1:nLS
                km1   = obj.getPreviousMaterialID(nLS,k);
                dt{k} = 0;
                for i = k:nLS
                    dt{k} = dt{k} -tfi{i}.*TD{i,km1} + tfi{km1}.*(tfiDer{i,k}.*TD{km1,i});
                end
            end
        end
    end

    methods (Access = private)
        function init(obj,desVar)
            obj.designVariable = desVar;
        end
    end

    methods (Static, Access = private)
        function km1 = getPreviousMaterialID(nLS,k)
            km1         = k-1;
            km1(km1==0) = nLS+1;
        end

        function xE = expandVariable(x,dJ) % Provisional until ndimf vector
            dJ = dJ.evaluate([0;0]);
            I  = size(x,1);
            J  = size(x,2);
            xE = cell(I,J);
            nE = ndims(dJ)-2;
            for i = 1:I
                for j = 1:J
                    xE{i,j} = Expand(x{i,j},nE);
                end
            end
        end
    end
end