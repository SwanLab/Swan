classdef MeshTypeComputer < handle
    
    properties (Access = private)
        nnode
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
                    switch obj.nnode
                        case 3
                            type = 'TRIANGLE';
                        case 4
                            type = 'QUAD';
                    end
                case 'Volume'
                    switch obj.nnode
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
            obj.nnode        =  cParams.nnode;
        end
        
    end
    
end

