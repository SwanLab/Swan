classdef RaviartThomasElement3D < handle
   
    properties (Access = private)
       
    end
    
    
    properties (Access = public)
        n_vertices
        vertices
        shapeFunctions
        normalVectors
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
            obj.computeNormalVectors()
            obj.computeShapeFunctions()
        end
        
        function computeVertices(obj)
            obj.n_vertices = 4;
            obj.vertices = [0,0,0;1,0,0;0,1,0;0,0,1];
        end
        
        function computeNormalVectors(obj)
            obj.normalVectors = [1,1,1;-1,0,0;0,-1,0;0,0,-1];
        end
        
        function F = integral_func(~,f,A,B,C)
            syms x y z a1 b1 a2 a3 t r real
            F = f;
            
            if A(1) == B(1) && A(1) == C(1)
                m = (A(3)-B(3))/(A(2)-B(2));
                n = A(3)-m*A(2);
                F = subs(F,x,A(1));
                F = subs(F,y,t);
                F = subs(F,z,m*t-n);
            elseif A(2) == B(2) && A(2) == C(2)
                m = (A(3)-B(3))/(A(1)-B(1));
                n = A(3)-m*A(1);
                F = subs(F,y,A(2));
                F = subs(F,x,t);
                F = subs(F,z,m*t-n);
            elseif A(3) == B(3) && A(3) == C(3)
                m = (A(2)-B(2))/(A(1)-B(1));
                n = A(2)-m*A(1);
                F = subs(F,y,m*t-n);
                F = subs(F,x,t);
                F = subs(F,z,A(3));
            else
                m = ((A(3)-B(3))/(A(2)-B(2))-(B(3)-C(3))/(B(2)-C(2)))/((A(1)-B(1))/(A(2)-B(2))-(B(1)-C(1))/(B(2)-C(2)));
                n = (A(3)-B(3))/(A(2)-B(2))-m*(A(1)-B(1))/(A(2)-B(2));
                l = A(3)-m*A(1)-n*A(2);
                F = subs(F,y,r);
                F = subs(F,x,t);
                F = subs(F,z,m*t+n*r+l);
            end
            
            f_int = F;
            F = int(f_int,t,0,1);
            
            F = int(F,r,0,1);
        end
        
        function computeShapeFunctions(obj)
            syms x y z a1 a2 a3 b1 real
            
            p = [a1+b1*x,a2+b1*y,a3+b1*z];
            for j = 1:obj.n_vertices
                pn(j) = dot(p,obj.normalVectors(j,:));
            end
            
            A(1) = obj.integral_func(pn(1),obj.vertices(2,:),obj.vertices(3,:),obj.vertices(4,:));
            A(2) = obj.integral_func(pn(2),obj.vertices(3,:),obj.vertices(1,:),obj.vertices(4,:));
            A(3) = obj.integral_func(pn(3),obj.vertices(1,:),obj.vertices(2,:),obj.vertices(4,:));
            A(4) = obj.integral_func(pn(4),obj.vertices(1,:),obj.vertices(2,:),obj.vertices(3,:));
            
            for i = 1:obj.n_vertices
                b = zeros(1,obj.n_vertices);
                b(i) = 1;
                
                eq = A == b;
                s = solve(eq,[a1 a2 a3 b1]);
                obj.shapeFunctions{i} = [s.a1+s.b1*x,s.a2+s.b1*y,s.a3+s.b1*z];
            end
        end
        
    end
    
end