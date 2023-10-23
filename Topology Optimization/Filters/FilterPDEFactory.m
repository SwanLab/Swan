classdef FilterPDEFactory < handle

    properties (Access = public)
        LHStype
        scale
    end

    methods (Access = public)

        function obj = FilterPDEFactory(cParams)
            metric = obj.evaluateMetric(cParams);
            if not(isfield(cParams,'boundaryType'))
                cParams.boundaryType = 'Neumann';
            end
            switch cParams.boundaryType
                case 'Neumann'
                    obj.scale = 'MACRO';
                    switch metric
                        case 'Isotropy'
                            obj.LHStype = 'StiffnessMass';
                        case 'Anisotropy'
                            obj.LHStype = 'AnisotropicStiffnessMass';
                    end
                case 'Robin'
                    obj.scale = 'MACRO';
                    switch metric
                        case 'Isotropy'
                            obj.LHStype = 'StiffnessMassBoundaryMass';
                        case 'Anisotropy'
                            obj.LHStype = 'AnisotropicStiffnessMassBoundaryMass';
                    end
                case 'Periodic'
                    obj.scale = 'MICRO';
                    switch metric
                        case 'Isotropy'
                            obj.LHStype = 'StiffnessMass';
                        case 'Anisotropy'
                            obj.LHStype = 'AnisotropicStiffnessMass';
                    end
            end
        end

    end

    methods (Access = private, Static)
        function metric = evaluateMetric(cParams)
            if isfield(cParams,'CAnisotropic')
                metric = 'Anisotropy';
            else
                metric = 'Isotropy';
            end
        end
    end

end