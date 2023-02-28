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
            for i=1:4
                subplot(2,2,i)
                hold on
                for j = 1:length(nodes)
                X(j,:) = obj.shapeFunctions{i}(nodes(j,1),nodes(j,2),nodes(j,3));
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
        
        function F = planeIntegral(~,func,pointA,pointB,pointC)
            syms x y z u v real
            paramx = (1-u-v)*pointA(1) + u*pointB(1) + v*pointC(1);
            paramy = (1-u-v)*pointA(2) + u*pointB(2) + v*pointC(2);
            paramz = (1-u-v)*pointA(3) + u*pointB(3) + v*pointC(3);
            
            func = subs(func,x,paramx);
            func = subs(func,y,paramy);
            func = subs(func,z,paramz);

            F = int(int(func,u,0,1),v,0,1);
        end
        
        function computeShapeFunctions(obj)
            syms x y z a1 a2 a3 b1 real
            
            p = [a1+b1*x,a2+b1*y,a3+b1*z];
            for j = 1:obj.n_vertices
                pn(j) = dot(p,obj.normalVectors(j,:));
            end
            
            A(1) = obj.planeIntegral(pn(1),obj.vertices(2,:),obj.vertices(3,:),obj.vertices(4,:));
            A(2) = obj.planeIntegral(pn(2),obj.vertices(3,:),obj.vertices(1,:),obj.vertices(4,:));
            A(3) = obj.planeIntegral(pn(3),obj.vertices(1,:),obj.vertices(2,:),obj.vertices(4,:));
            A(4) = obj.planeIntegral(pn(4),obj.vertices(1,:),obj.vertices(2,:),obj.vertices(3,:));
            
            for i = 1:obj.n_vertices
                b = zeros(1,obj.n_vertices);
                b(i) = 1;
                
                eq = A == b;
                s = solve(eq,[a1 a2 a3 b1]);
                obj.shapeFunctions{i} = matlabFunction([s.a1+s.b1*x,s.a2+s.b1*y,s.a3+s.b1*z]);
            end
        end
        
    end
    
end