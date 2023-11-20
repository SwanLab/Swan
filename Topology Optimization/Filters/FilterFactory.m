classdef FilterFactory < handle

    % - H-J: 'evaluate' method inside UnfittedBoundaryFunction !!
    % - try Clement, fVali=pi(Vi)

    methods (Access = public, Static)

        function filter = create(cParams)
            switch cParams.filterType
                case 'P1'
                    filter = FilterKernel(cParams);
                case 'AugmentedKernel'
                    filter = FilterKernelGeneral(cParams);
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
                    filter = FilterPDE(cParams);
            end
        end

    end

end