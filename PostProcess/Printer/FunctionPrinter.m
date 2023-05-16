classdef FunctionPrinter < handle
    
    methods (Static, Access = public)

        function obj = create(cParams)
            if isfield(cParams, 'type')
                switch cParams.type
                    case {'GiD', 'Gid', 'gid'}
                        obj = FunctionPrinter_GiD(cParams);
                    case {'Paraview'}
                        obj = FunctionPrinter_Paraview(cParams);
                end
            else
                % Default case
                obj = FunctionPrinter_Paraview(cParams);
            end
        end

    end

end

