classdef MaterialInterpolatorFactory < handle

    methods (Access = public, Static)

        function obj = create(cParams)
            switch cParams.interpolation
                case 'SIMPALL'
                    if ~isfield(cParams,'simpAllType')
                        cParams.simpAllType = 'EXPLICIT';
                    end
                    switch cParams.simpAllType
                        case 'EXPLICIT'
                            obj = SimpAllExplicitInterpolator(cParams);
                        case 'IMPLICIT'
                            switch cParams.dim
                                case '2D'
                                    obj = SimpAllInterpolationImplicit2D(cParams);
                                case '3D'
                                    obj = SimpAllInterpolationImplicit3D(cParams);
                                otherwise
                                    error('Invalid problem dimension.');
                            end
                        otherwise
                            error('Invalid SimpAll type.');
                    end
                case 'SIMP_Adaptative'
                    obj = SimpInterpolationAdaptative(cParams);
                case 'SIMP_P3'
                    obj = SimpInterpolationP3(cParams);
                case 'SIMPThermal'
                    obj = SIMPThermalInterpolation(cParams);
                case 'HomogenizedMicrostructure'
                    obj = HomogenizedMicrostructureInterpolator(cParams);
                case 'PhaseField'
                    if ~isfield(cParams,'subType')
                        cParams.subType = 'AT';
                    end
                    switch cParams.subType
                        case 'AT'
                            obj = PhaseFieldAmbrosioTortorelliDegradation(cParams);
                        case 'ATSplit'
                            obj = PhaseFieldAmbrosioTortorelliSplitDegradation(cParams);
                        case 'AT2linear'
                            obj = PhaseFieldAT2linearDegradation(cParams);
                        case 'Rational'
                            obj = PhaseFieldRationalDegradation(cParams);
                        case 'General'
                            obj = PhaseFieldGeneralDegradation(cParams);
                        case 'RationalWu'
                            obj = PhaseFieldRationalWuDegradation(cParams);
                    end
                otherwise
                    error('Invalid Material Interpolation method.');

            end

        end

    end
end