classdef MaterialFactory < handle

    methods (Access = public, Static)

        function material = create(cParams)

            switch cParams.type
                case 'ISOTROPIC'
                    switch cParams.ndim
                        case 2
                            material = Isotropic2dElasticMaterial(cParams);
                        case 3
                            material = Isotropic3dElasticMaterial(cParams);
                    end
                case 'ORTHOTROPIC'
                    material = Orthotropic2dElasticMaterial(cParams);
                    
                case 'HomogenizedMicrostructure'
                    material = HomogenizedMicrostructureInterpolator(cParams);
                    
                case 'DensityBased'
                    material = DensityBasedMaterial(cParams);
                
                case 'STOKES'
                    material = MaterialStokes(cParams);
                    
                case 'PhaseField'
                    switch cParams.PFtype
                        case 'Analytic'
                            material = MaterialPhaseFieldAnalytic(cParams);
                        case 'AnalyticSplit'
                            material = MaterialPhaseFieldAnalyticSplit(cParams);
                        case 'Homogenized'
                            material = MaterialPhaseFieldHomogenized(cParams);
                    end

                case 'ContinuumDamage'
                    material = MaterialContinuumDamage(cParams);
            end

        end

    end
    
end

