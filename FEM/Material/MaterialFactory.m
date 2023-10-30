classdef MaterialFactory < handle
    
    
    methods (Access = public, Static)
        
        function material = create(cParams)
            switch cParams.ptype
                case {'ELASTIC','DIFF-REACT'}
                    switch cParams.pdim
                        case '2D'
                            cParams.nstre = 3;    
                            material = Isotropic2dElasticMaterial(cParams);
                        case '3D'
                            cParams.nstre = 6;
                            material = Isotropic3dElasticMaterial(cParams);
                    end
                    
                case 'HYPERELASTIC'
                    error('Still not implemented.')
                    
                case 'THERMAL'
                    error('Still not implemented.')

                case 'Stokes'
                    material = Material_Stokes(cParams);
                    
                otherwise
                    error('Invalid ptype.')
            end
            
        end
        
        
    end
end