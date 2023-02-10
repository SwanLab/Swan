classdef CrouzeixRaviart1D < handle
   
    properties (Access = private)
       
    end
    
    
    properties (Access = public)
        n_vertices
        vertices
        n_nodes
        nodes
        edgeVectors
        midPoints
        shapeFunctions
        fig
    end
    
    
    methods (Access = public)
    
        function obj = CrouzeixRaviart1D()
            obj.init();
        end
        
        function plotShapeFunctions(obj)
            obj.fig = figure();
            syms x
            hold on
            for i=1:obj.n_vertices
                subplot(1,2,i)
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
       
        function init(obj)
            obj.computeVertices()
            obj.computeShapeFunctions()
        end
        
        function computeVertices(obj)
            obj.n_vertices = 2;
            obj.vertices = [0,1];
        end
        
        function computeShapeFunctions(obj)
            syms x
            for i = 1:obj.n_vertices
                A = [1 obj.vertices(2); 1 obj.vertices(1) ];
                X = zeros(obj.n_vertices,1);
                X(i) = 1;
                s = A\X;
                obj.shapeFunctions{i} = s(1)+s(2)*x;
            end
        end
        
    end
    
end