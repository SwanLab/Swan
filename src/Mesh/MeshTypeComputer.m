classdef MeshTypeComputer < handle
    
    methods (Static, Access = public)
        
        function type = compute(connec,g)
            nnodeElem = size(connec,2);
            switch g
                case 'Line'
                    type = 'LINE';
                case 'Surface'
                    switch nnodeElem
                        case 3
                            type = 'TRIANGLE';
                        case {4,9}
                            type = 'QUAD';
                        otherwise
                            type = 'SUPERELEMENT';
                    end
                case 'Volume'
                    switch nnodeElem
                        case 4
                            type = 'TETRAHEDRA';
                        case 8
                            type = 'HEXAHEDRA';
                    end
            end
        end
        
    end
    
end

