classdef MaterialFactory < handle
    
    
    methods (Access = public, Static)
        
        function material = create(cParams)
            switch cParams.ptype
                case {'ELASTIC','DIFF-REACT'}
                    switch cParams.pdim
                        case '2D'
                            material = Material_Elastic_ISO_2D(cParams);
                        case '3D'
                            material = Material_Elastic_ISO_3D(cParams);
                    end
                    
                case 'HYPERELASTIC'
                    switch cParams.pdim
                        case '2D'
                            cParams.connec = mesh.connec;
                            cParams.cartd  = geometry(1).cartd;
                            cParams.nnode  = geometry(1).interpolation.nnode;
                            cParams.coord  = mesh.coord;
                            material = Material_Hyperelastic_2D(cParams);
                    end
                    
                case 'THERMAL'
                    error('Still not implemented.')
                case 'Stokes'
                    material = Material_Stokes(nelem);
                otherwise
                    error('Invalid ptype.')
            end
            
        end
        
        
    end
end