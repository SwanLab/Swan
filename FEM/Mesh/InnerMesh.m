classdef InnerMesh < Mesh
    
    properties (Access = private)
        globalConnec
        backgroundCoord
    end
    
    methods (Access = public)
        
        function obj = InnerMesh(cParams)
            obj.init(cParams);
            obj.computeCoords();
            obj.computeConnec();
            obj.computeDescriptorParams();
            obj.createInterpolation();
            obj.computeElementCoordinates();            
            obj.unfittedType = 'SIMPLE';
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.globalConnec = cParams.globalConnec;
            obj.backgroundCoord = cParams.backgroundCoord;
        end
        
        function computeCoords(obj)
            coordElem = [];
            nnode = size(obj.globalConnec,2);
            for inode = 1:nnode
                nodes = obj.globalConnec(:,inode);
                coordElem = [coordElem; obj.backgroundCoord(nodes,:)];
            end
            subCoords = unique(coordElem,'rows','stable');
            obj.coord = subCoords;
        end
        
        function computeConnec(obj)
            connec = obj.globalConnec;
            coords = obj.backgroundCoord;
            subCoords = obj.coord;
            nnode = size(connec,2);
            for inode = 1:nnode
                coord = coords(connec(:,inode),:);
                I = obj.findIndexesComparingCoords(coord,subCoords);
                subConnec(:,inode) = I;
            end
            obj.connec = subConnec;
        end
        
    end
    
    methods (Access = private, Static)
        
        function I = findIndexesComparingCoords(A,B)
            I = zeros(1,size(A,1));
            for inode = 1:size(A,1)
                match = true(size(B,1),1);
                for idime = 1:size(A,2)
                    match = match & B(:,idime) == A(inode,idime);
                end
                I(inode) = find(match,1);
            end
        end
    end
    
end