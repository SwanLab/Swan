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
            for i=1:6
                subplot(2,3,i)
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
            
            baseShapeFunction = [a1+b3*y-b2*z,a2+b1*z-b3*x,a3+b2*x-b1*y];
            
            for j = 1:6
                tangentComponentBaseShapeFunction(j) = dot(baseShapeFunction,obj.edges.vect(j,:));
            end
            
            matrixLHS = obj.assemblyLHS(tangentComponentBaseShapeFunction);
            for i = 1:length(matrixLHS)
                vectorRHS = obj.assemblyRHS(i);
                coefShapeFunc = obj.computeShapeFunctionCoefficients(matrixLHS,vectorRHS);
                obj.shapeFunctions{i} = matlabFunction([coefShapeFunc.a1+coefShapeFunc.b3*y-coefShapeFunc.b2*z...
                    ,coefShapeFunc.a2+coefShapeFunc.b1*z-coefShapeFunc.b3*x,coefShapeFunc.a3+coefShapeFunc.b2*x...
                    -coefShapeFunc.b1*y],'Vars',[x y z]);
            end
        end
        
        function matrixLHS = assemblyLHS(obj,tangentComponentBaseShapeFunction)
            matrixLHS(1) = obj.lineIntegral(tangentComponentBaseShapeFunction(1),obj.vertices(2,:),obj.vertices(3,:));
            matrixLHS(2) = obj.lineIntegral(tangentComponentBaseShapeFunction(2),obj.vertices(2,:),obj.vertices(1,:));
            matrixLHS(3) = obj.lineIntegral(tangentComponentBaseShapeFunction(3),obj.vertices(1,:),obj.vertices(3,:));
            matrixLHS(4) = obj.lineIntegral(tangentComponentBaseShapeFunction(4),obj.vertices(1,:),obj.vertices(4,:));
            matrixLHS(5) = obj.lineIntegral(tangentComponentBaseShapeFunction(5),obj.vertices(2,:),obj.vertices(4,:));
            matrixLHS(6) = obj.lineIntegral(tangentComponentBaseShapeFunction(6),obj.vertices(3,:),obj.vertices(4,:));
        end
        
        function vectorRHS = assemblyRHS(obj,i)
            vectorRHS = zeros(1,length(obj.edges.vect));
            vectorRHS(i) = 1;
        end
        
        function s = computeShapeFunctionCoefficients(~,matrixLHS,vectorRHS)
            syms a1 a2 a3 b1 b2 b3 real
            eq = matrixLHS == vectorRHS;
            s = solve(eq,[a1 a2 a3 b1 b2 b3]);
        end
        
    end
    
end