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
            for i=1:3
                subplot(1,3,i)
                hold on
                for j = 1:length(nodes)
                    X(j,:) = obj.shapeFunctions{i}(nodes(j,1),nodes(j,2));
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
        
        function F = lineIntegral(~,func,pointA,pointB)
            syms x y a1 a2 b1 t real
            x1 = pointA(1); y1 = pointA(2);
            x2 = pointB(1); y2 = pointB(2);
            func = subs(func,x,x1 + t*(x2-x1));
            func = subs(func,y,y1 + t*(y2-y1));
            F = int(func,t,0,1);
        end
        
        function computeShapeFunctions(obj)
            syms x y a1 a2 b1 real
            
            baseShapeFunction = [a1+b1*x,a2+b1*y];
            for j = 1:obj.n_vertices
                normalComponentShapeFunction(j) = dot(baseShapeFunction,obj.normalVectors(j,:));
            end
            
            matrixLHS = obj.assemblyLHS(normalComponentShapeFunction);
            for i = 1:obj.n_vertices
                vectorRHS = obj.assemblyRHS(i);
                coefShapeFunc = solve(matrixLHS == vectorRHS,[a1 a2 b1]);
                obj.shapeFunctions{i} = matlabFunction([coefShapeFunc.a1+coefShapeFunc.b1*x,coefShapeFunc.a2+coefShapeFunc.b1*y]);
            end
        end
        
        function matrixLHS = assemblyLHS(obj,normalComponentShapeFunction)
            matrixLHS(1) = obj.lineIntegral(normalComponentShapeFunction(1),obj.vertices(2,:),obj.vertices(3,:));
            matrixLHS(2) = obj.lineIntegral(normalComponentShapeFunction(2),obj.vertices(3,:),obj.vertices(1,:));
            matrixLHS(3) = obj.lineIntegral(normalComponentShapeFunction(3),obj.vertices(1,:),obj.vertices(2,:));
        end
        
        function vectorRHS = assemblyRHS(obj,i)
            vectorRHS = zeros(1,obj.n_vertices);
            vectorRHS(i) = 1;
        end
        
    end
    
end

