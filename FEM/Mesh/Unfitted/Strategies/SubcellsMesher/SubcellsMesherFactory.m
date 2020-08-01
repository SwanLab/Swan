classdef SubcellsMesherFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            ndimIso = cParams.ndimIso;
            type    = cParams.type;
            switch ndimIso
                case 'Line'
                    obj = SubcellsMesher_1D(cParams);
                case 'Surface'
                    switch type
                        case 'INTERIOR'
                            obj  = SubcellsMesher_Interior(cParams);
                        case 'BOUNDARY'
                            obj = SubcellsMesher_Boundary_2D(cParams);
                    end
                case 'Volume'
                    switch type
                        case 'INTERIOR'
                            obj = SubcellsMesher_Interior(cParams);
                        case 'BOUNDARY'
                            obj = SubcellsMesher_Boundary_3D(cParams);
                    end
            end
        end
        
    end
end