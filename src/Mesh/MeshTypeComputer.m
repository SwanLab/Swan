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
                        case 20
                            type = 'HEXAHEDRA';  % Hexaedro cuadrático (20 nodos)
                        otherwise
                            % Para otros tipos de elementos volumétricos, usar HEXAHEDRA como default
                            % o TETRAHEDRA dependiendo del número de nodos
                            if nnodeElem > 8
                                type = 'HEXAHEDRA';  % Probablemente hexaedro de orden superior
                            else
                                type = 'TETRAHEDRA';  % Probablemente tetraedro de orden superior
                            end
                    end
            end
        end
        
    end
    
end

