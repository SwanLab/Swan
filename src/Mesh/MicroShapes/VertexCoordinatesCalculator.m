classdef VertexCoordinatesCalculator < handle
    
    properties (Access = private)
        latticeVectors
        nodes
    end
    
    properties (Access = public)
        vertCoord
    end
    
    methods (Access = public)
        
        function obj = VertexCoordinatesCalculator(cParams)
            obj.init(cParams);
            obj.computeVertex();
        end
        
        function computeVertex(obj)
            nVectors = size(obj.latticeVectors, 1);
            
            
            if nVectors == 2
                v1 = obj.latticeVectors(1,:);
                v2 = obj.latticeVectors(2,:);
                
                obj.vertCoord = [
                    0,    0;
                    v1;
                    v1 + v2;
                    v2
                ];
            else
                
                obj.obtainMasterVertex();
                obj.obtainSlaveVertex();
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.latticeVectors = cParams.latticeVectors;
            obj.nodes = cParams.nodes;
            obj.vertCoord = zeros(cParams.nodes.vert,2);
        end
        
        function obtainMasterVertex(obj)
            obj.vertCoord(1,:) = [0, 0];
            for iMaster = 1:obj.nodes.vert/2
                v = obj.latticeVectors(iMaster, :);
                obj.vertCoord(iMaster+1,:) = obj.vertCoord(iMaster,:) + v;
            end
        end
        
        function obtainSlaveVertex(obj)
            c0 = obj.vertCoord(obj.nodes.vert/2+1, :);
            obj.vertCoord(obj.nodes.vert/2+1,:) = c0;
            
            for iSlave = 1:obj.nodes.vert/2
                v = -obj.latticeVectors(iSlave, :);
                pos = c0 + v;
                
                if iSlave == obj.nodes.vert/2
                    if norm(obj.vertCoord(1,:) - pos) > 1e-10
                        warning('CRITICAL ERROR. Vertices computed wrongly');
                    end
                else
                    idx = obj.nodes.vert/2 + iSlave + 1;
                    obj.vertCoord(idx,:) = pos;
                    c0 = pos;
                end
            end
        end
        
    end

end