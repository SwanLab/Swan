classdef Simplicial1D < handle
   
    properties (Access = private)
    end
    
    
    properties (Access = public)
        vertices
        connectivities
        length
        tangentVectors
        normalVectors
    end
    
    
    methods (Access = public)
        
        function obj = Simplicial1D()
            obj.init();
        end
        
    end
    
    
    methods (Access = private)
        
        function init(obj)
            obj.computeVertices();
            obj.computeConnectivities();
            obj.computeLength();
            obj.computeTangentVectors();
            obj.computeNormalVectors();
            
        end
        
        function computeVertices(obj)
            obj.vertices = [0,1];
        end
        
        function computeConnectivities(obj)
            obj.connectivities = [1 2];
        end
        
        function computeLength(obj)
            obj.length = sqrt(sum((obj.vertices(2)-obj.vertices(1)).^2));
        end
        
        function computeTangentVectors(obj)
            obj.tangentVectors = (obj.vertices(2)-obj.vertices(1))/obj.length;
        end
        
        function computeNormalVectors(obj)
            obj.normalVectors = obj.tangentVectors * [1 0;0 -1];
        end
        
    end
    
end