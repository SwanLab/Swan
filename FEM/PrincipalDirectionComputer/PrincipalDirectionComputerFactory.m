classdef PrincipalDirectionComputerFactory < handle
    
    methods (Access = public, Static)
        
        function pdc = create(cParams)
            
            switch cParams.type
                case '2D'
                    pdc = PrincipalDirectionComputerIn2D(cParams);
                case '3D'
                    pdc = PrincipalDirectionComputerIn3D(cParams);
            end
        end
        
    end
    
end
