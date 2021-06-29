classdef MeshPlotterFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            ndimIso = cParams.ndimIso;
            type    = cParams.type;
            switch ndimIso
                case 1
                    obj = MeshPlotter_1D;
                case 2
                    switch type
                        case 'INTERIOR'
                            obj = MeshPlotter_Interior_2D;
                        case 'BOUNDARY'
                            obj = MeshPlotter_Boundary_2D;
                    end
                case 3
                    switch type
                        case 'INTERIOR'
                            obj = MeshPlotter_Null;
                        case 'BOUNDARY'
                            obj = MeshPlotter_Boundary_3D;
                    end
            end
            
        end
        
    end
    
    
end