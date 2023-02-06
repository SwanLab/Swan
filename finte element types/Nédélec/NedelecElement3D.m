classdef NedelecElement3D < handle
   
    properties (Access = private)
       
    end
    
    
    properties (Access = public)
        n_vertices
        vertices
        edges
        midPoints
        shapeFunctions
        fig
    end
    
    
    methods (Access = public)
    
        function obj = NedelecElement3D()
            obj.init();
        end
        
        function plotShapeFunctions(obj)
            obj.fig = figure();
            nodes = [0,0,0;0,1/3,0;0,2/3,0;0,1,0;1/3,0,0;1/3,1/3,0;1/3,2/3,0;2/3,0,0;2/3,1/3,0;2/3,1/3,0;1,0,0;
                     0,0,1/3;0,1/3,1/3;0,2/3,1/3;1/3,0,1/3;1/3,1/3,1/3;2/3,0,1/3;0,0,2/3;0,1/3,2/3;1/3,0,2/3;0,0,1];
            syms x y z
            for i=1:6
                subplot(2,3,i)
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
            obj.computeEdges()
            obj.computeMidPoints()
            obj.computeShapeFunctions()
        end
        
        function computeVertices(obj)
            obj.n_vertices = 4;
            obj.vertices = [0,0,0;0,1,0;1,0,0;0,0,1];
        end
        
        function computeEdges(obj)
            obj.edges.vect = [sqrt(2)/2,-sqrt(2)/2,0;-1,0,0;0,1,0;0,0,1;0,-sqrt(2)/2,sqrt(2)/2;-sqrt(2)/2,0,sqrt(2)/2];
            obj.edges.measure = [sqrt(2),1,1,1,sqrt(2),sqrt(2)];
            obj.edges.j = [4,5,6,1,2,3];
        end
        
        function computeMidPoints(obj)
            obj.midPoints = [0.5,0.5,0;0.5,0,0;0,0.5,0;0,0,0.5;0,0.5,0.5;0.5,0,0.5];
        end
        
        function computeShapeFunctions(obj)
            syms x y z
            for i = 1:6
                obj.shapeFunctions{i} = cross([x y z]-obj.midPoints(obj.edges.j(i),:),obj.edges.vect(obj.edges.j(i),:))/...
                    (obj.edges.measure(i)*dot(obj.edges.vect(i,:),cross(obj.midPoints(i,:)-obj.midPoints(obj.edges.j(i),:),obj.edges.vect(obj.edges.j(i),:))));
            end
        end
        
    end
    
end