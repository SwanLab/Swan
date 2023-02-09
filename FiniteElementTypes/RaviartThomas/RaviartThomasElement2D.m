classdef RaviartThomasElement2D < handle
   
    properties (Access = private)
       
    end
    
    
    properties (Access = public)
        n_vertices
        vertices
        measure
        shapeFunctions
        fig
    end
    
    
    methods (Access = public)
    
        function obj = RaviartThomasElement2D()
            obj.init();
        end
        
        function plotShapeFunctions(obj)
            obj.fig = figure();
            nodes = [0,0;0,1/3;0,2/3;0,1;1/3,0;1/3,1/3;1/3,2/3;2/3,0;2/3,1/3;2/3,1/3;1,0];
            syms x y
            for i=1:3
                subplot(1,3,i)
                hold on
                for j = 1:length(nodes)
                    X(j,:) = subs(obj.shapeFunctions{i},[x y],nodes(j,:));
                end
                quiver(nodes(:,1),nodes(:,2),double(X(:,1)),double(X(:,2)))
                plot([0 1],[0 0],'k')
                plot([0 0],[1 0],'k')
                plot([0 1],[1 0],'k')
                title("i:"+i)
                xlabel('x')
                ylabel('y')
                grid on
                hold off
            end
        end
    
    end
    
    
    methods (Access = private)
       
        function init(obj)
            obj.computeVertices()
            obj.computeMeasure()
            obj.computeShapeFunctions()
        end
        
        function computeVertices(obj)
            obj.n_vertices = 3;
            obj.vertices = [0,0;0,1;1,0];
        end
        
        function computeMeasure(obj)
            obj.measure = (obj.vertices(2,2)*obj.vertices(3,1))/2;
        end
        
        function computeShapeFunctions(obj)
            for i = 1:obj.n_vertices
                syms x y
                obj.shapeFunctions{i} = 1/(2*obj.measure)*([x,y]-obj.vertices(i,:));
            end
        end
        
    end
    
end