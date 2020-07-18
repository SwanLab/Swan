classdef MeshTypeComputer < handle
    
    properties (Access = private)
       ndim
       nnode
       kFace
    end        
    
    methods (Access = public)
        
        function obj = MeshTypeComputer(cParams)
            obj.init(cParams)
        end
        
        function type = compute(obj)            
            switch obj.ndim
                case 1
                    type = 'LINE';
                case 2
                    switch obj.kFace
                        case -1
                            type = 'LINE';
                        case 0
                            switch obj.nnode
                                case 2
                                    type = 'LINE';
                                case 3
                                    type = 'TRIANGLE';
                                case 4
                                    type = 'QUAD';
                            end
                    end
                case 3
                    switch obj.kFace
                        case -2
                            type = 'LINE';
                        case -1
                            switch obj.nnode
                                case 3
                                    type = 'TRIANGLE';
                                case 4
                                    type = 'QUAD';
                            end
                        case 0
                            switch obj.nnode                             
                                case 4
                                    type = 'TETRAHEDRA';
                                case 8
                                    type = 'HEXAHEDRA';
                            end
                    end
            end
            if obj.nnode == 0
                type = 'EMPTY MESH';
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams) 
            obj.ndim  = cParams.ndim;
            obj.nnode = cParams.nnode;
            obj.kFace = cParams.kFace;
        end
        
    end
    
end

