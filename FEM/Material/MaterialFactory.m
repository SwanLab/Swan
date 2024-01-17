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
                case 'ANISOTROPIC'
                            material = AnisotropicGiven(cParams);
                    


            end

        end

    end
    
end

