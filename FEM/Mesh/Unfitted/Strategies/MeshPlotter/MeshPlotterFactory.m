classdef MeshPlotterFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            ndimIso      = cParams.ndimIso;
            unfittedType = cParams.unfittedType;
            switch ndimIso
                case 1
                    obj = MeshPlotter_1D;
                case 2
                    switch unfittedType
                        case 'INTERIOR'
                            obj = MeshPlotter_Interior_2D;
                        case 'BOUNDARY'
                            obj = MeshPlotter_Boundary_2D;
                    end
                case 3
                    switch unfittedType
                        case 'INTERIOR'
                            obj = MeshPlotter_Null;
                        case 'BOUNDARY'
                            obj = MeshPlotter_Boundary_3D;
                    end
            end
            
        end
        
    end
    
    
end