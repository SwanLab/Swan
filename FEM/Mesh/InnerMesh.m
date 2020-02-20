classdef InnerMesh < Mesh
    
    properties (Access = private)
        globalConnec
        backgroundCoord
        isInBoundary
    end
    
    methods (Access = public)
        
        function obj = InnerMesh(cParams)
            obj.init(cParams);
            obj.computeCoords();
            obj.computeConnec();
            obj.computeDescriptorParams();
            obj.createInterpolation();
            obj.computeElementCoordinates();
        end
        
        function add2plot(obj,ax)
            patch(ax,'vertices',obj.coord,'faces',obj.connec,...
                'edgecolor',[0.5 0 0], 'edgealpha',0.5,'edgelighting','flat',...
                'facecolor',[1 0 0],'facelighting','flat')
            axis('equal');
        end
        
    end
    
    methods (Access = protected)
        
        function computeEmbeddingDim(obj)
            if obj.isInBoundary
                obj.embeddedDim = obj.ndim - 1;
            else
                obj.embeddedDim = obj.ndim;
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.globalConnec = cParams.globalConnec;
            obj.backgroundCoord = cParams.backgroundCoord;
            obj.unfittedType = 'SIMPLE';
            obj.isInBoundary = cParams.isInBoundary;
            if isfield(cParams,'type')
                obj.type = cParams.type;
            else
                obj.type = 'INTERIOR';
            end
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