classdef InnerMesh < Mesh
    
    properties (GetAccess = public, SetAccess = private)
        backgroundMesh        
    end
    
    properties (Access = private)
        globalConnec
        all2unique
        unique2all
        uniqueNodes
    end
    
    methods (Access = public)
        
        function obj = InnerMesh(cParams)
            obj.init(cParams);
            obj.computeUniqueNodes();
            obj.computeCoords();
            obj.computeConnec();
            obj.computeDescriptorParams();
            obj.createInterpolation();
            obj.computeElementCoordinates();
        end
        

        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.globalConnec   = cParams.globalConnec;
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.isInBoundary   = cParams.isInBoundary;
            obj.type = 'INTERIOR';
        end
        
        function computeUniqueNodes(obj)
            allNodes = obj.globalConnec(:);
            [uNodes,ind,ind2] = unique(allNodes,'rows','stable');
            obj.all2unique  = ind;    
            obj.unique2all  = ind2;   
            obj.uniqueNodes = uNodes;
        end
        
        function computeCoords(obj)
            uNodes       = obj.uniqueNodes;
            allCoords    = obj.backgroundMesh.coord;
            uniqueCoords = allCoords(uNodes,:);
            obj.coord    = uniqueCoords;
        end
        
        function computeConnec(obj)
            nnode = size(obj.globalConnec,2);
            nCell = size(obj.globalConnec,1);             
            obj.connec = reshape(obj.unique2all,nCell,nnode);
        end          
        
    end
    
end