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
                    
                case 'HomogenizedMicrostructure'
                    material = HomogenizedMicrostructureInterpolator(cParams);
                    
                case 'DensityBased'
                    material = DensityBasedMaterial(cParams);
                
                case 'STOKES'
                    material = Material_Stokes(cParams);
                    
                case 'PhaseFieldAnalytic'
                    material = MaterialPhaseField(cParams);
                    
                case 'PhaseFieldHomog'
                    material = HomogenizedPhaseField(cParams);
            end

        end

    end
    
end

