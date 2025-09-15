classdef UnfittedMeshFunction < handle

    methods (Static, Access = public)

        function obj = create(cParams)
            f = cParams.fun;
            switch class(f)
                case 'LagrangianFunction'
                    obj = UnfittedMeshLagrangian(cParams);
                    obj.compute(f);
                case 'AnalyticalFunction'
                    obj = UnfittedMeshAnalytical(cParams);
                    obj.compute(f);
            end
        end

    end

end