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
            obj.computeShapeFunctions()
        end
        
        function computeVertices(obj)
            obj.n_vertices = 3;
            obj.vertices = [0,0;1,0;0,1];
        end
        
        function computeEdges(obj)
            obj.edges.vect = [-1,1;0,-1;1,0];
            obj.edges.measure = [sqrt(2),1,1];
        end
        
        function F = lineIntegral(~,func,pointA,pointB)
            syms x y t real
            x1 = pointA(1); y1 = pointA(2);
            x2 = pointB(1); y2 = pointB(2);
            func = subs(func,x,x1 + t*(x2-x1));
            func = subs(func,y,y1 + t*(y2-y1));
            F = int(func,t,0,1);
        end
        
        function computeShapeFunctions(obj)
            syms x y a1 a2 b1 real
            
            p = [a1-b1*y,a2+b1*x];
            for j = 1:obj.n_vertices
                pn(j) = dot(p,obj.edges.vect(j,:));
            end
            
            A(1) = obj.lineIntegral(pn(1),obj.vertices(2,:),obj.vertices(3,:));
            A(2) = obj.lineIntegral(pn(2),obj.vertices(3,:),obj.vertices(1,:));
            A(3) = obj.lineIntegral(pn(3),obj.vertices(1,:),obj.vertices(2,:));
            
            for i = 1:obj.n_vertices
                b = zeros(1,obj.n_vertices);
                b(i) = 1;
                
                eq = A == b;
                s = solve(eq,[a1 a2 b1]);
                obj.shapeFunctions{i} = [s.a1-s.b1*y,s.a2+s.b1*x];
            end
        end
        
    end
    
end