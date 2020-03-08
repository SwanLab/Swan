classdef InnerMesh < Mesh
    
    properties (Access = private)
        globalConnec
        backgroundCoord
        isInBoundary
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
        
        function computeUniqueNodes(obj)
            allNodes = obj.globalConnec(:);
            [uNodes,ind,ind2] = unique(allNodes,'rows','stable');
            obj.all2unique  = ind;    
            obj.unique2all  = ind2;   
            obj.uniqueNodes = uNodes;
        end
        
        function computeCoords(obj)
            uNodes       = obj.uniqueNodes;
            allCoords    = obj.backgroundCoord;
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