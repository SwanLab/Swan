classdef Simplicial3D < handle
   
    properties (Access = private)
    end
    
    
    properties (Access = public)
        vertices
        connectivities
        edgesLength
        facesArea
        volume
        tangentVectors
        normalVectors
    end
    
    
    methods (Access = public)
        
        function obj = Simplicial3D()
            obj.init();
        end
        
    end
    
    
    methods (Access = private)
        
        function init(obj)
            obj.computeVertices();
            obj.computeConnectivities();
            obj.computeEdgesLength();
            obj.computeTangentVectors();
            obj.computeNormalVectors();
            obj.computeFacesArea();
            obj.computeVolume();
            
        end
        
        function computeVertices(obj)
            obj.vertices = [0,0,0;1,0,0;0,1,0;0,0,1];
        end
        
        function computeConnectivities(obj)
            obj.connectivities = [1 2 3 4];
        end
        
        function computeEdgesLength(obj)
            len(1) = sqrt(sum((obj.vertices(2,:)-obj.vertices(3,:)).^2));
            len(2) = sqrt(sum((obj.vertices(1,:)-obj.vertices(3,:)).^2));
            len(3) = sqrt(sum((obj.vertices(1,:)-obj.vertices(2,:)).^2));
            len(4) = sqrt(sum((obj.vertices(4,:)-obj.vertices(1,:)).^2));
            len(5) = sqrt(sum((obj.vertices(4,:)-obj.vertices(2,:)).^2));
            len(6) = sqrt(sum((obj.vertices(4,:)-obj.vertices(3,:)).^2));
            
            obj.edgesLength = len;
        end
        
        function computeFacesArea(obj)
            area(1) = obj.edgesLength(2)*obj.edgesLength(3)*0.5;
            area(2) = obj.edgesLength(2)*obj.edgesLength(4)*0.5;
            area(3) = obj.edgesLength(3)*obj.edgesLength(4)*0.5;
            area(4) = 0.5*dot(obj.edgesLength(1).*obj.tangentVectors(1,:),obj.edgesLength(5).*obj.tangentVectors(5,:));
            
            obj.facesArea = area;
        end
        
        function computeVolume(obj)
            obj.volume = obj.edgesLength(2)*obj.edgesLength(3)*obj.edgesLength(4)*0.5;
        end
        
        function computeTangentVectors(obj)
            tang(1,:) = (obj.vertices(3,:)-obj.vertices(2,:))/obj.edgesLength(1);
            tang(2,:) = (obj.vertices(1,:)-obj.vertices(3,:))/obj.edgesLength(2);
            tang(3,:) = (obj.vertices(2,:)-obj.vertices(1,:))/obj.edgesLength(3);
            tang(4,:) = (obj.vertices(4,:)-obj.vertices(1,:))/obj.edgesLength(4);
            tang(5,:) = (obj.vertices(4,:)-obj.vertices(2,:))/obj.edgesLength(5);
            tang(6,:) = (obj.vertices(4,:)-obj.vertices(3,:))/obj.edgesLength(6);
            
            obj.tangentVectors = tang;
        end
        
        function computeNormalVectors(obj)
            norm(1,:) = cross(obj.tangentVectors(5,:),obj.tangentVectors(6,:));
            norm(2,:) = cross(obj.tangentVectors(2,:),obj.tangentVectors(4,:));
            norm(3,:) = cross(obj.tangentVectors(3,:),obj.tangentVectors(4,:));
            norm(4,:) = cross(obj.tangentVectors(3,:),obj.tangentVectors(2,:));
            
            norm(1,:) = norm(1,:)/(sqrt(sum(norm(1,:).^2)));
            norm(2,:) = norm(2,:)/(sqrt(sum(norm(2,:).^2)));
            norm(3,:) = norm(3,:)/(sqrt(sum(norm(3,:).^2)));
            norm(4,:) = norm(4,:)/(sqrt(sum(norm(4,:).^2)));
            
            obj.normalVectors = norm;
        end
        
    end
    
end