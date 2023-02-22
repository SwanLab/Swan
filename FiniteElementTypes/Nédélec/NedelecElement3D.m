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
            obj.computeShapeFunctions()
        end
        
        function computeVertices(obj)
            obj.n_vertices = 4;
            obj.vertices = [0,0,0;1,0,0;0,1,0;0,0,1];
        end
        
        function computeEdges(obj)
            obj.edges.vect = [1,-1,0;1,0,0;0,1,0;0,0,1;-1,0,1;0,1,-1];
            obj.edges.measure = [sqrt(2),1,1,1,sqrt(2),sqrt(2)];
            obj.edges.j = [4,5,6,1,2,3];
        end
        
        function F = integral_func(~,f,A,B)
            syms x y z a1 a2 a3 b1 b2 b3 t real
            F = f;
            
            if A(1) == B(1)
                F = subs(F,x,A(1));
                if A(2) == B(2)
                    F = subs(F,y,A(2));
                    F = subs(F,z,t);
                elseif A(3) == B(3)
                    F = subs(F,y,t);
                    F = subs(F,z,A(3));
                else
                    m = (A(3)-B(3))/(A(2)-B(2));
                    n = A(3)-m*A(2);
                    F = subs(F,y,t);
                    F = subs(F,z,m*t+n);
                end
            elseif A(2) == B(2)
                F = subs(F,y,A(2));
                if A(3) == B(3)
                    F = subs(F,z,A(3));
                    F = subs(F,x,t);
                elseif A(1) ~= B(1)
                    m = (A(3)-B(3))/(A(1)-B(1));
                    n = A(3)-m*A(1);
                    F = subs(F,x,t);
                    F = subs(F,z,m*t+n);
                end
            elseif A(3) == B(3)
                m = (A(2)-B(2))/(A(1)-B(1));
                n = A(2)-m*A(1);
                F = subs(F,y,m*t+n);
                F = subs(F,x,t);
                F = subs(F,z,A(3));
            else
                m1 = (A(2)-B(2))/(A(1)-B(1));
                n1 = A(2)-m*A(1);
                m2 = (A(3)-B(3))/(A(1)-B(1));
                n2 = A(3)-m*A(1);
                F = subs(F,x,t);
                F = subs(F,y,m1*t+n1);
                F = subs(F,z,m2*t+n2);
            end
            
            f_int = F;
            F = int(f_int,t,0,1);
        end
        
        function F = lineIntegral(~,func,pointA,pointB)
            syms x y z t real
            x1 = pointA(1); y1 = pointA(2); z1 = pointA(3);
            x2 = pointB(1); y2 = pointB(2); z2 = pointB(3);
            func = subs(func,x,x1 + t*(x2-x1));
            func = subs(func,y,y1 + t*(y2-y1));
            func = subs(func,z,z1 + t*(z2-z1));
            F = int(func,t,0,1);
        end
        
        function computeShapeFunctions(obj)
            syms x y z a1 a2 a3 b1 b2 b3 real
            
            xx = cross([x y z],[b1 b2 b3]);
            p = [a1+b3*y-b2*z,a2+b1*z-b3*x,a3+b2*x-b1*y];
            for j = 1:6
                pn(j) = dot(p,obj.edges.vect(j,:));
            end
            
            A(1) = obj.lineIntegral(pn(1),obj.vertices(2,:),obj.vertices(3,:));
            A(2) = obj.lineIntegral(pn(2),obj.vertices(2,:),obj.vertices(1,:));
            A(3) = obj.lineIntegral(pn(3),obj.vertices(1,:),obj.vertices(3,:));
            A(4) = obj.lineIntegral(pn(4),obj.vertices(1,:),obj.vertices(4,:));
            A(5) = obj.lineIntegral(pn(5),obj.vertices(2,:),obj.vertices(4,:));
            A(6) = obj.lineIntegral(pn(6),obj.vertices(3,:),obj.vertices(4,:));
            
            for i = 1:length(A)
                b = zeros(1,length(A));
                b(i) = 1;
                
                eq = A == b;
                s = solve(eq,[a1 a2 a3 b1 b2 b3]);
                obj.shapeFunctions{i} = [s.a1+s.b3*y-s.b2*z,s.a2+s.b1*z-s.b3*x,s.a3+s.b2*x-s.b1*y];
            end
        end
        
    end
    
end