classdef VertexCoordinatesCalculator < handle
    
    properties (Access = private)
        c
        theta
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
            obj.obtainMasterVertex();
            obj.obtainSlaveVertex();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.c = cParams.c;
            obj.theta = cParams.theta;
            obj.nodes = cParams.nodes;
            obj.vertCoord = zeros(cParams.nodes.vert,2);
        end
        
        function obtainMasterVertex(obj)
            c0 = [0,0];
            for iMaster = 1:obj.nodes.vert/2
                leng = obj.c(iMaster);
                angle = obj.theta(iMaster);
                pos = NodeCoordinatesComputer.computeThePosition(c0,leng,angle);
                obj.vertCoord(iMaster+1,:) = obj.vertCoord(iMaster+1,:)+pos;
                c0 = pos;
            end
        end
        
        function obtainSlaveVertex(obj)
            c0 = obj.vertCoord(obj.nodes.vert/2+1,:);
            for iSlave = 1:obj.nodes.vert/2
                leng = obj.c(iSlave);
                angle = obj.theta(iSlave)+180;
                pos = NodeCoordinatesComputer.computeThePosition(c0,leng,angle);
                if iSlave == obj.nodes.vert/2
                    if obj.vertCoord(1,:) ~= pos
                        cprintf('red','CRYTICAL ERROR. Vertices computed wrongly \n');
                    end
                else
                    iMaster = obj.nodes.vert/2;
                    obj.vertCoord(iMaster+iSlave+1,:) = obj.vertCoord(iMaster+iSlave+1,:)+pos;
                    c0 = pos;
                end
            end
        end
        
    end

end