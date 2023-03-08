classdef RaviartThomasElement2D < handle
   
    properties (Access = private)
       
    end
    
    
    properties (Access = public)
        ndofs
        vertices
        normalVectors
        shapeFunctions
        fig
        simplicial
    end
    
    
    methods (Access = public)
    
        function obj = RaviartThomasElement2D()
            obj.init();
        end
        
        function plotShapeFunctions(obj)
            obj.fig = figure();
            nodes = [0,0;0,1/3;0,2/3;0,1;1/3,0;1/3,1/3;1/3,2/3;2/3,0;2/3,1/3;2/3,1/3;1,0];
            for i=1:3
                subplot(1,3,i);
                hold on
                for j = 1:length(nodes)
                    X(j,:) = obj.shapeFunctions{i}(nodes(j,1),nodes(j,2));
                end
                quiver(nodes(:,1),nodes(:,2),double(X(:,1)),double(X(:,2)));
                plot([0 1],[0 0],'k'); plot([0 0],[1 0],'k'); plot([0 1],[1 0],'k');
                title("i:"+i);
                xlabel('x'); ylabel('y');
                grid on
                hold off
            end
        end
    
    end
    
    
    methods (Access = private)
       
        function init(obj)
            obj.simplicial = Simplicial2D();
            obj.computeVertices();
            obj.computeNdof();
            obj.computeNormalVectors();
            obj.computeShapeFunctions();
        end
        
        function computeVertices(obj)
            obj.vertices = obj.simplicial.vertices;
        end
        
        function computeNdof(obj)
            obj.ndofs = length(obj.vertices);
        end
        
        function computeNormalVectors(obj)
            obj.normalVectors = obj.simplicial.normalVectors;
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
            syms x y
            LHS = obj.applyLinearForms();
            for s = 1:obj.ndofs
                RHS = obj.computeLinearFormsValues(s);
                c = obj.computeShapeFunctionCoefficients(LHS,RHS);
                obj.shapeFunctions{s} = matlabFunction([c.a1+c.b1*x,c.a2+c.b1*y]);
            end
        end
        
        function c = computeShapeFunctionCoefficients(~,LHS,RHS)
            syms a1 a2 b1
            c = solve(LHS == RHS,[a1 a2 b1]);
        end
        
        function LHS = applyLinearForms(obj,normalComponentShapeFunction)
            syms x y a1 a2 b1 real
            
            baseShapeFunction = [a1+b1*x,a2+b1*y];
            for j = 1:obj.ndofs
                normalComponentShapeFunction(j) = dot(baseShapeFunction,obj.normalVectors(j,:));
            end
            
            v = obj.vertices;
            LHS(1) = obj.lineIntegral(normalComponentShapeFunction(1),v(2,:),v(3,:));
            LHS(2) = obj.lineIntegral(normalComponentShapeFunction(2),v(3,:),v(1,:));
            LHS(3) = obj.lineIntegral(normalComponentShapeFunction(3),v(1,:),v(2,:));
        end
        
        function RHS = computeLinearFormsValues(obj,s)
            RHS = zeros(1,obj.ndofs);
            RHS(s) = 1;
        end
        
    end
    
end

