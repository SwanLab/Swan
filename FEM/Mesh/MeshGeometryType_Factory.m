classdef MeshGeometryType_Factory < handle
    
    methods (Access = public, Static)
        
        function geometryType = getGeometryType(ndim,nnode,embeddedDim)
            switch ndim
                case 1
                    geometryType = 'LINE';
                case 2
                    switch embeddedDim
                        case 1
                            geometryType = 'LINE';
                        case 2
                            switch nnode
                                case 2
                                    geometryType = 'LINE';
                                case 3
                                    geometryType = 'TRIANGLE';
                                case 4
                                    geometryType = 'QUAD';
                            end
                    end
                case 3
                    switch embeddedDim
                        case 1
                            geometryType = 'LINE';
                        case 2
                            switch nnode
                                case 3
                                    geometryType = 'TRIANGLE';
                                case 4
                                    geometryType = 'QUAD';
                            end
                        case 3
                            switch nnode                             
                                case 4
                                    geometryType = 'TETRAHEDRA';
                                case 8
                                    geometryType = 'HEXAHEDRA';
                            end
                    end
            end
            if nnode == 0
                geometryType = 'EMPTY MESH';
            end
        end
        
    end
    
end

