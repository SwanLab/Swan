classdef FilterFactory < handle

    % Refactoring ToDo:
    % General:
    % - GeometricalFunction
    % - delete: computeRHSInBoundaries, FilterPDEUnfitted
    % - create: UnfittedFunction, BoundaryUnfittedFunction
    % - 'evaluate' method inside RHSUnfitted
    % - Test P0Function in unfitted not available

    methods (Access = public, Static)

        function filter = create(cParams)
            switch cParams.filterType
                case 'P1'
                    switch cParams.designVarType
                        case {'Continuous','Density','MicroParams'}
                            filter = FilterP1(cParams);
                        case 'LevelSet'
                            filter = FilterP1Unfitted(cParams);
                    end
                case 'PDE'
                    if not(isfield(cParams,'boundaryType'))
                        cParams.boundaryType = 'Neumann';
                    end
                    if not(isfield(cParams,'metric'))
                        cParams.metric       = 'Isotropy';
                    end
                    switch cParams.boundaryType
                        case {'Neumann','Periodic'}
                            switch cParams.metric
                                case 'Isotropy'
                                    cParams.LHStype = 'StiffnessMass';
                                case 'Anisotropy'
                                    cParams.LHStype = 'AnisotropicStiffnessMass';
                            end
                        case 'Robin'
                            switch cParams.metric
                                case 'Isotropy'
                                    cParams.LHStype = 'StiffnessMassBoundaryMass';
                                case 'Anisotropy'
                                    cParams.LHStype = 'AnisotropicStiffnessMassBoundaryMass';
                            end
                    end
                    switch cParams.designVarType
                        case {'Continuous','Density','MicroParams'}
                            filter = FilterPDE(cParams);
                        case 'LevelSet'
                            filter = FilterPDEUnfitted(cParams);
                    end
            end
        end

    end

end