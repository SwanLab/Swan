classdef SubcellsMesherFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            ndimIso      = cParams.ndimIso;
            unfittedType = cParams.unfittedType;
            switch ndimIso
                case 1
                    obj = SubcellsMesher_1D;
                case 2
                    switch unfittedType
                        case 'INTERIOR'
                            obj  = SubcellsMesher_Interior;
                        case 'BOUNDARY'
                            obj = SubcellsMesher_Boundary_2D;
                    end
                case 3
                    switch unfittedType
                        case 'INTERIOR'
                            obj = SubcellsMesher_Interior;
                        case 'BOUNDARY'
                            obj = SubcellsMesher_Boundary_3D;
                    end
            end
        end
        
    end
end