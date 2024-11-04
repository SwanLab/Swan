classdef MultimaterialGradientComputer < handle

    properties (Access = private)
        mesh
        designVariable
    end

    methods (Access = public)
        function obj = MultimaterialGradientComputer(cParams)
            obj.init(cParams);
        end

        function dt = compute(obj,TD)
            x       = obj.designVariable;
            ls2     = x.levelSets{1,2};
            ls3     = x.levelSets{1,3};
            tfi  = x.obtainDomainFunction();
            chi2 = obj.computeExactCharFun(ls2);
            chi3 = obj.computeExactCharFun(ls3);

            dt{1} = - tfi{1}.*TD{1,end} - tfi{2}.*TD{2,end} - tfi{3}.*TD{3,end} ...
                + tfi{4}.*( (1-chi2).*TD{end,1} + (1-chi3).*chi2.*TD{end,2} + chi2.*chi3.*TD{end,3} );

            dt{2} = - tfi{2}.*TD{2,1} - tfi{3}.*TD{3,1} + tfi{1}.*( (1-chi3).*TD{1,2} + chi3.*TD{1,3} );

            dt{3} = tfi{2}.*TD{2,3} - tfi{3}.*TD{3,2};
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh           = cParams.mesh;
            obj.designVariable = cParams.designVariable;
        end

        function f = computeUnitAnalyticalFunction(obj)
            s.fHandle = @(x) ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            f         = AnalyticalFunction(s);
        end

        function chi = computeExactCharFun(obj,ls)
            chiLS   = ls.getUnfittedMesh();
            s.uMesh = chiLS;
            s.fun   = obj.computeUnitAnalyticalFunction();
            uFun    = UnfittedFunction(s);
            chi     = uFun.project('P0');
        end

    end

end