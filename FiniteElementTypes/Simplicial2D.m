classdef Simplicial2D < handle
   
    properties (Access = private)
    end
    
    
    properties (Access = public)
        vertices
        connectivities
        edgesLength
        area
        tangentVectors
        normalVectors
    end
    
    
    methods (Access = public)
        
        function obj = Simplicial2D()
            obj.init();
        end
        
    end
    
    
    methods (Access = private)
        
        function init(obj)
            obj.computeVertices();
            obj.computeConnectivities();
            obj.computeEdgesLength;
            obj.computeAera();
            obj.computeTangentVectors();
            obj.computeNormalVectors();
            
        end
        
        function computeVertices(obj)
            obj.vertices = [0,0;0,1;1,0];
        end
        
        function computeConnectivities(obj)
            obj.connectivities = [1 2 3];
        end
        
        function computeEdgesLength(obj)
            len(1) = sqrt(sum((obj.vertices(2,:)-obj.vertices(3,:)).^2));
            len(2) = sqrt(sum((obj.vertices(1,:)-obj.vertices(3,:)).^2));
            len(3) = sqrt(sum((obj.vertices(1,:)-obj.vertices(2,:)).^2));
            
            obj.edgesLength = len;
        end
        
        function computeAera(obj)
            obj.area = obj.edgesLength(2)*obj.edgesLength(3)*0.5;
        end
        
        function computeTangentVectors(obj)
            tang(1,:) = (obj.vertices(3,:)-obj.vertices(2,:))/obj.edgesLength(1);
            tang(2,:) = (obj.vertices(1,:)-obj.vertices(3,:))/obj.edgesLength(2);
            tang(3,:) = (obj.vertices(2,:)-obj.vertices(1,:))/obj.edgesLength(3);
            
            obj.tangentVectors = tang;
        end
        
        function computeNormalVectors(obj)
            obj.normalVectors = obj.tangentVectors * [1 0;0 -1];
        end
        
    end
    
end