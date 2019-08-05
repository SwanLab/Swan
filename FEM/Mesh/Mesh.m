classdef Mesh < AbstractMesh & matlab.mixin.Copyable
    
    properties (GetAccess = public, SetAccess = protected)
%         coord
%         connec
        ndim
        nnode
        npnod
%         nelem
    end
    
    properties (Access = public)
       % unfittedType
        meshBackground
    end
    
    methods (Access = public)
        
        function obj = create(obj,coord,connec)
            obj.coord  = coord;
            obj.connec = connec;
            obj.computeDescriptorParams();
            obj.unfittedType = 'SIMPLE';
        end
        
        function objClone = clone(obj)
            objClone = copy(obj);
        end
        
        function S = computeMeanCellSize(obj)
            x1 = obj.coord(obj.connec(:,1));
            x2 = obj.coord(obj.connec(:,2));
            x3 = obj.coord(obj.connec(:,3));
            
            x1x2 = abs(x2-x1);
            x2x3 = abs(x3-x2);
            x1x3 = abs(x1-x3);
            hs = max([x1x2,x2x3,x1x3]');
            
            S = mean(hs);
        end
        
        function L = computeCharacteristicLength(obj)
            xmin = min(obj.coord);
            xmax = max(obj.coord);
            L = norm(xmax-xmin)/2;
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
            obj.computeGeometryType();
        end
        
    end
    
    methods (Access = private)
        
        function computeGeometryType(obj)
            factory = MeshGeometryType_Factory();
            obj.geometryType = factory.getGeometryType(obj.ndim,obj.nnode);
        end
        
    end
    
end

