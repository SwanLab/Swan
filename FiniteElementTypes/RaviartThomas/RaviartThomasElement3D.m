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
            
            baseShapeFunction = [a1+b1*x,a2+b1*y,a3+b1*z];
            for j = 1:obj.n_vertices
                normalComponentShapeFunction(j) = dot(baseShapeFunction,obj.normalVectors(j,:));
            end
            
            matrixLHS = obj.assemblyLHS(normalComponentShapeFunction);
            for i = 1:obj.n_vertices
                vectorRHS = obj.assemblyRHS(i);
                coefShapeFunc = solve(matrixLHS == vectorRHS,[a1 a2 a3 b1]);
                obj.shapeFunctions{i} = matlabFunction([coefShapeFunc.a1+coefShapeFunc.b1*x,coefShapeFunc.a2+coefShapeFunc.b1*y,coefShapeFunc.a3+coefShapeFunc.b1*z]);
            end
        end
        
        function matrixLHS = assemblyLHS(obj,normalComponentShapeFunction)
            matrixLHS(1) = obj.planeIntegral(normalComponentShapeFunction(1),obj.vertices(2,:),obj.vertices(3,:),obj.vertices(4,:));
            matrixLHS(2) = obj.planeIntegral(normalComponentShapeFunction(2),obj.vertices(3,:),obj.vertices(1,:),obj.vertices(4,:));
            matrixLHS(3) = obj.planeIntegral(normalComponentShapeFunction(3),obj.vertices(1,:),obj.vertices(2,:),obj.vertices(4,:));
            matrixLHS(4) = obj.planeIntegral(normalComponentShapeFunction(4),obj.vertices(1,:),obj.vertices(2,:),obj.vertices(3,:));
        end
        
        function vectorRHS = assemblyRHS(obj,i)
            vectorRHS = zeros(1,obj.n_vertices);
            vectorRHS(i) = 1;
        end
        
    end
    
end