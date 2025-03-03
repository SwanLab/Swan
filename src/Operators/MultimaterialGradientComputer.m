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
            nLS    = length(tfi)-1;
            dt     = cell(nLS,1);
            for k = 1:nLS
                km1   = obj.getPreviousMaterialID(nLS,k);
                dt{k} = 0;
                for i = k:nLS
                    dt{k} = dt{k} -DDP(tfi{i},TD{i,km1}) + DDP(tfi{km1},DDP(tfiDer{i,k},TD{km1,i}));
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
    end
end