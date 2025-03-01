classdef MeshTypeComputer < handle
    
    properties (Access = private)
        nnodeElem
        geometryType
    end
    
    methods (Access = public)
        
        function obj = MeshTypeComputer(cParams)
            obj.init(cParams)
        end
        
        function type = compute(obj)
            switch obj.geometryType
                case 'Line'
                    type = 'LINE';
                case 'Surface'
                    switch obj.nnodeElem
                        case 3
                            type = 'TRIANGLE';
                        case 4
                            type = 'QUAD';
                    end
                case 'Volume'
                    switch obj.nnodeElem
                        case 4
                            type = 'TETRAHEDRA';
                        case 8
                            type = 'HEXAHEDRA';
                    end
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.geometryType = cParams.geometryType;
            obj.nnodeElem    = cParams.nnodeElem;
        end
        
    end
    
end

