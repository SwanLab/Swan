classdef SplitterInConnectedComponents < handle
    
    properties (Access = private)
        edges
        adjancyMatrix
        graphV
        connectedComp
    end
    
    properties (Access = private)
        faces
    end
    
    methods (Access = public)
        
        function obj = SplitterInConnectedComponents(cParams)
            obj.init(cParams)
        end
        
        function c = split(obj)
            obj.computeEdges();
            obj.computeAdjacencyMatrix();
            obj.computeGraph();
            obj.computeConnectedComponents();
            c = obj.connectedComp;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.faces = cParams.faces;
        end
        
        function computeEdges(obj)
            s.type = 'TRIANGLE';
            s.nodesByElem = obj.faces;
            edge = EdgesConnectivitiesComputer(s);
            edge.compute();
            obj.edges = edge.nodesInEdges;
        end
        
        function [A] = computeAdjacencyMatrix(obj)
            nodeA = obj.edges(:,1);
            nodeB = obj.edges(:,2);
            A = sparse(nodeA,nodeB,1);
            obj.adjancyMatrix = A+A';
        end
        
        function computeGraph(obj)
            A = obj.adjancyMatrix;
            g = graph(A,'omitselfloops');
            obj.graphV = g;
        end
        
        function computeConnectedComponents(obj)
            g = obj.graphV;
            comp = conncomp(g);
            obj.connectedComp = comp;
        end
    
    end

end