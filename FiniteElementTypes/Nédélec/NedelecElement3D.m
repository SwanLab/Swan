classdef NedelecElement3D < handle
   
    properties (Access = private)
       
    end
    
    
    properties (Access = public)
        vertices
        ndofs
        edges
        shapeFunctions
        fig
        simplicial
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
                quiver3(nodes(:,1),nodes(:,2),nodes(:,3),double(X(:,1)),double(X(:,2)),double(X(:,3)));
                plot3([0 1],[0 0],[0 0],'k'); plot3([0 0],[1 0],[0 0],'k'); plot3([0 1],[1 0],[0 0],'k');
                plot3([0 0],[0 0],[0 1],'k'); plot3([0 0],[1 0],[0 1],'k'); plot3([1 0],[0 0],[0 1],'k');
                xlabel('x'); ylabel('y'); zlabel('z')
                title("i:"+string(i-1))
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
            obj.computeEdges();
            obj.computeShapeFunctions();
        end
        
        function computeVertices(obj)
            obj.vertices = obj.simplicial.vertices;
        end
        
        function computeNdof(obj)
            obj.ndofs = length(obj.simplicial.edgesLength);
        end
        
        function computeEdges(obj)
            obj.edges.vect = obj.simplicial.tangentVectors;
            obj.edges.measure = obj.simplicial.edgesLength;
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
            syms x y z
            LHS = obj.applyLinearForms();
            for s = 1:length(LHS)
                RHS = obj.computeLinearFormValues(s);
                c = obj.computeShapeFunctionCoefficients(LHS,RHS);
                obj.shapeFunctions{s} = matlabFunction([c.a1+c.b3*y-c.b2*z,c.a2+c.b1*z-c.b3*x,c.a3+c.b2*x-c.b1*y],'Vars',[x y z]);
            end
        end
        
        function LHS = applyLinearForms(obj)
            syms x y z a1 a2 a3 b1 b2 b3 real
            
            baseShapeFunction = [a1+b3*y-b2*z,a2+b1*z-b3*x,a3+b2*x-b1*y];
            
            for j = 1:obj.ndofs
                tangCompBaseShapeFunction(j) = dot(baseShapeFunction,obj.edges.vect(j,:));
            end
            
            v = obj.vertices;
            LHS(1) = obj.lineIntegral(tangCompBaseShapeFunction(1),v(2,:),v(3,:));
            LHS(2) = obj.lineIntegral(tangCompBaseShapeFunction(2),v(3,:),v(1,:));
            LHS(3) = obj.lineIntegral(tangCompBaseShapeFunction(3),v(2,:),v(1,:));
            LHS(4) = obj.lineIntegral(tangCompBaseShapeFunction(4),v(1,:),v(4,:));
            LHS(5) = obj.lineIntegral(tangCompBaseShapeFunction(5),v(2,:),v(4,:));
            LHS(6) = obj.lineIntegral(tangCompBaseShapeFunction(6),v(3,:),v(4,:));
        end
        
        function RHS = computeLinearFormValues(obj,s)
            RHS = zeros(1,length(obj.edges.vect));
            RHS(s) = obj.edges.measure(s);
        end
        
        function c = computeShapeFunctionCoefficients(~,LHS,RHS)
            syms a1 a2 a3 b1 b2 b3 real
            eq = LHS == RHS;
            c = solve(eq,[a1 a2 a3 b1 b2 b3]);
        end
        
    end
    
end