classdef Mesh < AbstractMesh & matlab.mixin.Copyable
    
    properties (GetAccess = public, SetAccess = protected)
        nnode
        npnod

        embeddedDim
    end
    
    properties (Access = public)
        meshBackground
    end    
       
    properties (Access = protected)
       type 
    end
    
    methods (Access = public)
        
        function obj = create(obj,cParams)
            obj.init(cParams);
            obj.computeDescriptorParams();
            obj.createInterpolation();
            obj.computeElementCoordinates();
        end
        
        function objClone = clone(obj)
            objClone = copy(obj);
        end
        
        function plot(obj)
            figure;
            patch('vertices',obj.coord,'faces',obj.connec,...
                'edgecolor',[0.5 0 0], 'edgealpha',0.5,'edgelighting','flat',...
                'facecolor',[1 0 0],'facelighting','flat')
            axis('equal');
        end
        
        function S = computeMeanCellSize(obj)
            x1(:,1) = obj.coord(obj.connec(:,1),1);
            x1(:,2) = obj.coord(obj.connec(:,1),2);
            x2(:,1) = obj.coord(obj.connec(:,2),1);
            x2(:,2) = obj.coord(obj.connec(:,2),2);
            x3(:,1) = obj.coord(obj.connec(:,3),1);
            x3(:,2) = obj.coord(obj.connec(:,3),2);            
            x1x2 = (x2-x1);
            x2x3 = (x3-x2);
            x1x3 = (x1-x3);
            n12 = sqrt(x1x2(:,1).^2 + x1x2(:,2).^2);
            n23 = sqrt(x2x3(:,1).^2 + x2x3(:,2).^2);
            n13 = sqrt(x1x3(:,1).^2 + x1x3(:,2).^2);
            hs = max([n12,n23,n13],[],2);
            S = max(hs);
        end
        
        function L = computeCharacteristicLength(obj)
            xmin = min(obj.coord);
            xmax = max(obj.coord);
            L = norm(xmax-xmin);%/2;
        end
        
        function changeCoordinates(obj,newCoords)
            obj.coord = newCoords;
        end
        
        function setCoord(obj,newCoord)
            obj.coord = newCoord;
        end
        
        function setConnec(obj,newConnec)
            obj.connec = newConnec;
        end
        
    end
    
    methods (Access = protected)
        
        function computeDescriptorParams(obj)
            obj.npnod = size(obj.coord ,1);
            obj.ndim  = size(obj.coord, 2);
            obj.nelem = size(obj.connec,1);
            obj.nnode = size(obj.connec,2);
            obj.computeEmbeddingDim();
            obj.computeGeometryType();
        end
        
        function computeEmbeddingDim(obj)
            switch obj.type
                case 'BOUNDARY'
                    obj.embeddedDim = obj.ndim - 1;
                case {'INTERIOR','COMPOSITE'}
                    obj.embeddedDim = obj.ndim;
                otherwise
                    error('EmbeddedDim not defined')
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.coord  = cParams.coord;
            obj.connec = cParams.connec;
            if isfield(cParams,'type')
                obj.type = cParams.type;
            else
                obj.type = 'INTERIOR';
            end
            obj.unfittedType = 'SIMPLE';
            
        end
        
        function computeGeometryType(obj)
            factory = MeshGeometryType_Factory();
            obj.geometryType = factory.getGeometryType(obj.ndim,obj.nnode,obj.embeddedDim);
        end
        
    end
    
end

