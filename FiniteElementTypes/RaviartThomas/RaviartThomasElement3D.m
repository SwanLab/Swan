classdef RaviartThomasElement3D < handle
   
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
    
        function obj = RaviartThomasElement3D()
            obj.init();
        end
        
        function plotShapeFunctions(obj)
            obj.fig = figure();
            nodes = [0,0,0;0,1/3,0;0,2/3,0;0,1,0;1/3,0,0;1/3,1/3,0;1/3,2/3,0;2/3,0,0;2/3,1/3,0;2/3,1/3,0;1,0,0;
                     0,0,1/3;0,1/3,1/3;0,2/3,1/3;1/3,0,1/3;1/3,1/3,1/3;2/3,0,1/3;0,0,2/3;0,1/3,2/3;1/3,0,2/3;0,0,1];
            syms x y z
            for i=1:4
                subplot(2,2,i)
                hold on
                for j = 1:length(nodes)
                X(j,:) = subs(obj.shapeFunctions{i},[x y z],nodes(j,:));
                end
                quiver3(nodes(:,1),nodes(:,2),nodes(:,3),double(X(:,1)),double(X(:,2)),double(X(:,3)))
                plot3([0 1],[0 0],[0 0],'k')
                plot3([0 0],[1 0],[0 0],'k')
                plot3([0 1],[1 0],[0 0],'k')
                plot3([0 0],[0 0],[0 1],'k')
                plot3([0 0],[1 0],[0 1],'k')
                plot3([1 0],[0 0],[0 1],'k')
                xlabel('x')
                ylabel('y')
                zlabel('z')
                title("i:"+i)
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
            obj.n_vertices = 4;
            obj.vertices = [0,0,0;0,1,0;1,0,0;0,0,1];
        end
        
        function computeMeasure(obj)
            obj.measure = (obj.vertices(2,2)*obj.vertices(3,1)*obj.vertices(4,3))/2;
        end
        
        function computeShapeFunctions(obj)
            for i = 1:obj.n_vertices
                syms x y z
                obj.shapeFunctions{i} = 1/(2*obj.measure)*([x,y,z]-obj.vertices(i,:));
            end
        end
        
    end
    
end