classdef FunctionPrinter < handle
    
    methods (Static, Access = public)

        function obj = create(cParams)
            if ~isfield(cParams,'type')
                run('UserVariables.m')
                cParams.type = default_export_software;
            end

            if isfield(cParams, 'type')
                switch cParams.type
                    case {'GiD', 'Gid', 'gid', 'GID'}
                        obj = FunctionPrinterGiD(cParams);
                    case {'Paraview', 'paraview', 'pvw', 'pv'}
                        obj = FunctionPrinterParaview(cParams);
                end
            else
                % Default case
                obj = FunctionPrinter_Paraview(cParams);
            end
        end

    end

end

