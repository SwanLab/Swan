classdef ConnectedVertexesComputer < handle
    
    properties (Access = private)
       edges 
    end
    
    methods (Access = public)
        
        function obj = ConnectedVertexesComputer(cParams)
            obj.init(cParams)            
        end
        
        function v = compute(obj,vertex)
            otherVertexWhenA  = obj.computeVertexBWhenVertexIsA(vertex);   
            otherVertexWhenB  = obj.computeVertexAWhenVertexIsB(vertex);   
            v = [otherVertexWhenA;otherVertexWhenB];                        
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.edges = cParams.edges;
        end
                    
        function otherVertex = computeVertexBWhenVertexIsA(obj,currentVertex)
            localCurrentVertex = 1;
            localOtherVertex   = 2;
            lC = localCurrentVertex;
            lO = localOtherVertex;
            cV = currentVertex;
            oV = obj.computeOtherVertexOfEdge(cV,lC,lO); 
            otherVertex = oV;                       
        end
        
        function otherVertex = computeVertexAWhenVertexIsB(obj,currentVertex)
            localCurrentVertex = 2;
            localOtherVertex   = 1;
            lC = localCurrentVertex;
            lO = localOtherVertex;
            cV = currentVertex;
            oV = obj.computeOtherVertexOfEdge(cV,lC,lO); 
            otherVertex = oV;
        end        
                
        function otherVertex = computeOtherVertexOfEdge(obj,currentVertex,localCurrentVertex,localOtherVertex)
            vertexInEdges        = obj.edges.nodesInEdges;            
            hasEdgeCurrentVertex = vertexInEdges(:,localCurrentVertex) == currentVertex;
            otherVertex          = vertexInEdges(hasEdgeCurrentVertex,localOtherVertex);            
        end             
        
    end
    
end