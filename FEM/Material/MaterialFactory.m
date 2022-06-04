classdef MaterialFactory < handle
    
    
    methods (Access = public, Static)
        
        function material = create(cParams)
            switch cParams.ptype
                case {'ELASTIC','DIFF-REACT','EIGENMODES'}
                    switch cParams.pdim
                        case '2D'
                            cParams.nstre = 3;    
                            material = Isotropic2dElasticMaterial(cParams);
                        case '3D'
                            cParams.nstre = 6;
                            material = Isotropic3dElasticMaterial(cParams);
                    end
                    
                case 'HYPERELASTIC'
                    switch cParams.pdim
                        case '2D'
                            cParams.connec = cParams.mesh.connec;
                            cParams.dNdx   = cParams.geometry.dNdx;
                            cParams.nnode  = cParams.mesh.nnodeElem;
                            cParams.coord  = cParams.mesh.coord;
                            material = Isotropic2dHyperElasticMaterial(cParams);
                    end
                    
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