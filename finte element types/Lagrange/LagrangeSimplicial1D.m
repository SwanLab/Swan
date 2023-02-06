classdef LagrangeSimplicial1D < handle
   
    properties (Access = private)
        k
    end
    
    
    properties (Access = public)
        n_vertices
        vertices
        normalVectors
        n_nodes
        nodes
        barycentricCoords
        shapeFunctions
        fig
    end
    
    
    methods (Access = public)
    
        function obj = LagrangeSimplicial1D(k)
            obj.init(k);
        end
    
        function plotShapeFunctions(obj)
            obj.fig = figure();
            hold on
            for i = 1:obj.k+1
                subplot(1,obj.k+1,i)
                fplot(obj.shapeFunctions{i},[0 1])
                xlabel('x')
                ylabel('y')
                title("i:"+i)
                grid on
            end
            hold off
        end
        
    end
    
    
    methods (Access = private)
       
        function init(obj,k)
            obj.k = k;
            obj.computeVertices()
            obj.computeNormalVectors()
            obj.computeNodes()
            obj.computeBarycentricCoords()
            obj.computeShapeFunctions()
        end
        
        function computeVertices(obj)
            obj.n_vertices = 2;
            obj.vertices = [0,1];
        end
        
        function computeNormalVectors(obj)
            obj.normalVectors = [1,0];
        end
        
        function computeNodes(obj)
            obj.n_nodes = nchoosek(1+obj.k,obj.k);
            
            if obj.k == 1
                obj.nodes.coord = [0,1];
                obj.nodes.index = [1,2];
            elseif obj.k == 2
                obj.nodes.coord = [0,0.5,1];
                obj.nodes.index = [1,0;1,2;2,0];
            elseif obj.k == 3
                obj.nodes.coord = [0,1/3,2/3,1];
                obj.nodes.index = [1,0;1,2;2,1;2,0];
            end
        end
        
        function computeBarycentricCoords(obj)
            syms x
            for i = 1:obj.n_vertices
                if i~=obj.n_vertices
                    obj.barycentricCoords{i} = symfun(1-(x-obj.vertices(i))/(obj.vertices(i+1)-obj.vertices(i)),x);
                else
                    obj.barycentricCoords{i} = symfun(1-(x-obj.vertices(i))/(obj.vertices(1)-obj.vertices(i)),x);
                end
            end
        end
        
        function computeShapeFunctions(obj)
            if obj.k==1
                obj.shapeFunctions = obj.barycentricCoords;
                
            elseif obj.k==2
                obj.shapeFunctions = cell(obj.n_nodes,1);
                for i = 1:obj.n_nodes
                    if obj.nodes.index(i,2) == 0
                        obj.shapeFunctions{i} = obj.barycentricCoords{obj.nodes.index(i,1)}*(2*obj.barycentricCoords{obj.nodes.index(i,1)}-1);
                    else
                        obj.shapeFunctions{i} = 4*obj.barycentricCoords{obj.nodes.index(i,1)}*obj.barycentricCoords{obj.nodes.index(i,2)};
                    end
                end
                
            elseif obj.k==3
                obj.shapeFunctions = cell(obj.n_nodes,1);
                for i = 1:obj.n_nodes
                    if obj.nodes.index(i,2) == 0
                        func = obj.barycentricCoords{obj.nodes.index(i,1)};
                        obj.shapeFunctions{i} = 0.5*func*(3*func-1)*(3*func-2);
                    else 
                        func_i = obj.barycentricCoords{obj.nodes.index(i,1)};
                        func_j = obj.barycentricCoords{obj.nodes.index(i,2)};
                        obj.shapeFunctions{i} = 9/2*func_i*(3*func_i-1)*func_j;
                    end
                end
            end
        end
        
    end
    
end