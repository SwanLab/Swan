classdef RaviartThomasElement3D < handle
   
    properties (Access = private)
       
    end
    
    
    properties (Access = public)
        vertices
        shapeFunctions
        normalVectors
        fig
        simplicial
        ndofs
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
                quiver3(nodes(:,1),nodes(:,2),nodes(:,3),double(X(:,1)),double(X(:,2)),double(X(:,3)));
                plot3([0 1],[0 0],[0 0],'k'); plot3([0 0],[1 0],[0 0],'k'); plot3([0 1],[1 0],[0 0],'k');
                plot3([0 0],[0 0],[0 1],'k'); plot3([0 0],[1 0],[0 1],'k'); plot3([1 0],[0 0],[0 1],'k');
                xlabel('x'); ylabel('y'); zlabel('z');
                title("i: "+string(i-1))
                grid on
                hold off
            end
        end
    
    end
    
    
    methods (Access = private)
       
        function init(obj)
            obj.simplicial = Simplicial3D();
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
            for j = 1:obj.ndofs
                normalComponentShapeFunction(j) = dot(baseShapeFunction,obj.normalVectors(j,:));
            end
            
            LHS = obj.applyLinearForm(normalComponentShapeFunction);
            for s = 1:obj.ndofs
                RHS = obj.computeLinearFormValues(s);
                c = solve(LHS == RHS,[a1 a2 a3 b1]);
                obj.shapeFunctions{s} = matlabFunction([c.a1+c.b1*x,c.a2+c.b1*y,c.a3+c.b1*z]);
            end
        end
        
        function LHS = applyLinearForm(obj,func)
            v = obj.vertices;
            LHS(1) = obj.planeIntegral(func(1),v(2,:),v(3,:),v(4,:));
            LHS(2) = obj.planeIntegral(func(2),v(3,:),v(1,:),v(4,:));
            LHS(3) = obj.planeIntegral(func(3),v(1,:),v(2,:),v(4,:));
            LHS(4) = obj.planeIntegral(func(4),v(1,:),v(2,:),v(3,:));
        end
        
        function RHS = computeLinearFormValues(obj,s)
            RHS = zeros(1,obj.ndofs);
            RHS(s) = obj.simplicial.facesArea(s);
        end
        
    end
    
end