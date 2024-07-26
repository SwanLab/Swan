classdef FilterFactory < handle

    methods (Access = public, Static)

        function filter = create(cParams)
            switch cParams.filterType
                case 'P1'
                    filter = FilterKernel(cParams);
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
                        case 'DirichletProjection'
                            cParams.LHStype = 'MassBoundaryMass';
                    end
                    filter = FilterPDE(cParams);
                case 'LUMP'
                    filter = FilterLump(cParams);
                case 'Filter&Project'
                    filter = FilterAndProject(cParams);
            end
        end

    end

end