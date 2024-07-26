classdef CharacteristicFunction < handle

    methods (Static, Access = public)
        function obj = create(uMesh)
            m             = uMesh.backgroundMesh;
            fHandle       = @(x) ones(size(x(1,:,:)));
            s.ndimf       = 1;
            s.fHandle     = fHandle;
            s.mesh        = m;
            aFun          = AnalyticalFunction(s);
            cParams.fun   = aFun;
            cParams.uMesh = uMesh;
            obj           = UnfittedFunction(cParams);
        end

        function obj = createAtBoundary(uMesh)
            m             = uMesh.backgroundMesh;
            fHandle       = @(x) ones(size(x(1,:,:)));
            s.ndimf       = 1;
            s.fHandle     = fHandle;
            s.mesh        = m;
            aFun          = AnalyticalFunction(s);
            cParams.fun   = aFun;
            cParams.uMesh = uMesh;
            obj           = UnfittedBoundaryFunction(cParams);
        end
    end

end