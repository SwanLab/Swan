classdef PrincipalDirectionComputerFactory < handle
    
    methods (Access = public, Static)
        
        function pdc = create(cParams)
            
            switch cParams.type
                case 2
                    pdc = PrincipalDirectionComputerIn2D(cParams);
                case 3
                    pdc = PrincipalDirectionComputerIn3D(cParams);
            end
        end
        
    end
    
end
