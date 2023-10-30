classdef CharacteristicFunction < handle

    methods (Static, Access = public)
        function obj = create(cParams)
            m           = cParams.levelSet.mesh;
            fHandle     = @(x) 1;
            s.ndimf     = 1;
            s.fHandle   = fHandle;
            s.mesh      = m;
            aFun        = AnalyticalFunction(s);
            cParams.fun = aFun;
            obj         = UnfittedFunction(cParams);
        end
    end

end