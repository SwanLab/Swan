classdef NedelecElement2D < handle
   
    properties (Access = private)
       
    end
    
    
    properties (Access = public)
        n_vertices
        vertices
        edges
        measure
        shapeFunctions
        fig
    end
    
    
    methods (Access = public)
    
        function obj = NedelecElement2D()
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
            obj.computeEdges()
            obj.computeMeasure()
            obj.computeShapeFunctions()
        end
        
        function computeVertices(obj)
            obj.n_vertices = 3;
            obj.vertices = [0,0;0,1;1,0];
        end
        
        function computeEdges(obj)
            obj.edges.vect = [sqrt(2)/2,-sqrt(2)/2;-1,0;0,1];
            obj.edges.measure = [sqrt(2),1,1];
        end
        
        function computeMeasure(obj)
            obj.measure = (obj.vertices(2,2)*obj.vertices(3,1))/2;
        end
        
        function [b] = Rmap(obj,a)
            b(1) = -a(2);
            b(2) = a(1);
        end
        
        function computeShapeFunctions(obj)
                syms x y
                obj.shapeFunctions{1} = obj.Rmap([x y]-obj.vertices(1,:))/(obj.edges.measure(1)*...
                    dot(obj.edges.vect(1,:),obj.Rmap((obj.vertices(2,:)+obj.vertices(3,:))/2-obj.vertices(1,:))));
                obj.shapeFunctions{2} = obj.Rmap([x y]-obj.vertices(2,:))/(obj.edges.measure(2)*...
                    dot(obj.edges.vect(2,:),obj.Rmap((obj.vertices(1,:)+obj.vertices(3,:))/2-obj.vertices(2,:))));
                obj.shapeFunctions{3} = obj.Rmap([x y]-obj.vertices(3,:))/(obj.edges.measure(3)*...
                    dot(obj.edges.vect(3,:),obj.Rmap((obj.vertices(2,:)+obj.vertices(1,:))/2-obj.vertices(3,:))));

        end
        
    end
    
end