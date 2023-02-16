classdef RaviartThomasElement2D < handle
   
    properties (Access = private)
       
    end
    
    
    properties (Access = public)
        n_vertices
        vertices
        measure
        normalVectors
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
            obj.computeNormalVectors()
            obj.computeShapeFunctions()
        end
        
        function computeVertices(obj)
            obj.n_vertices = 3;
            obj.vertices = [0,0;1,0;0,1];
        end
        
        function computeNormalVectors(obj)
            obj.normalVectors = [1,1;-1,0;0,-1];
        end
        
        function F = integral_func(~,f,A,B)
            syms x y a1 b1 a2 b2 real
            F = f;
            if A(1) ~= B(1)
                f_int = F;
                b = max(A(1),B(1));
                a = min(A(1),B(1));
                F = int(f_int,x,a,b);
            else
                f_int = F;
                F = subs(f_int,x,A(1));
            end
            if A(2) ~= B(2)
                f_int = F;
                b = max(A(2),B(2));
                a = min(A(2),B(2));
                F = int(f_int,y,a,b);
            else
                f_int = F;
                F = subs(f_int,y,A(2));
            end
        end
        
        function computeShapeFunctions(obj)
            syms x y a1 a2 b1 real
            
            p = [a1+b1*x,a2+b1*y];
            for j = 1:obj.n_vertices
                pn(j) = dot(p,obj.normalVectors(j,:));
            end
            
            A(1) = obj.integral_func(pn(1),obj.vertices(2,:),obj.vertices(3,:));
            A(2) = obj.integral_func(pn(2),obj.vertices(3,:),obj.vertices(1,:));
            A(3) = obj.integral_func(pn(3),obj.vertices(1,:),obj.vertices(2,:));
            
            for i = 1:obj.n_vertices
                b = zeros(1,obj.n_vertices);
                b(i) = 1;
                
                eq = A == b;
                s = solve(eq,[a1 a2 b1]);
                obj.shapeFunctions{i} = [s.a1+s.b1*x,s.a2+s.b1*y];
            end
        end
        
    end
    
end

