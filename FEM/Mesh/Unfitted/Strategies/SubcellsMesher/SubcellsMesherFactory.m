classdef SubcellsMesherFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            ndimIso = cParams.ndimIso;
            type    = cParams.type;
            switch ndimIso
                case 1
                    obj = SubcellsMesher_1D(cParams);
                case 2
                    switch type
                        case 'INTERIOR'
                            obj  = SubcellsMesher_Interior(cParams);
                        case 'BOUNDARY'
                            obj = SubcellsMesher_Boundary_2D(cParams);
                    end
                case 3
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