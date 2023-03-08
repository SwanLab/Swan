classdef CrouzeixRaviart2D < handle
   
    properties (Access = private)
       
    end
    
    
    properties (Access = public)
        vertices
        ndofs
        shapeFunctions
        fig
        simplicial
    end
    
    
    methods (Access = public)
    
        function obj = CrouzeixRaviart2D()
            obj.init();
        end
        
        function plotShapeFunctions(obj)
            set(groot,'defaulttextinterpreter','latex');
            set(groot,'defaultLegendInterpreter','latex');
            set(groot,'defaultAxesTickLabelInterpreter','latex');
            
            s.coord = obj.vertices;
            s.connec = [1 2 3];
            m = Mesh(s);
            
%             obj.fig = figure();
            for i=1:obj.ndofs
%                 subplot(obj.polinomialOrder+1,obj.polinomialOrder+1,i);
                figure()
                m.plot();
                trisurf(m.connec,m.coord(:,1),m.coord(:,2),obj.shapeFunctions{i}(m.coord(:,1),m.coord(:,2)));
                
                xlim([0 1]); ylim([0 1]); zlim([-1.5 1.5]);
                xlabel('x'); ylabel('y'); zlabel('z');
                title("i:"+string(i-1));
                grid on
            end
        end
    
    end
    
    
    methods (Access = private)
       
        function init(obj)
            obj.simplicial = Simplicial2D();
            obj.computeVertices();
            obj.computeNdof();
            obj.computeShapeFunctions();
        end
        
        function computeVertices(obj)
            obj.vertices = obj.simplicial.vertices;
        end
        
        function computeNdof(obj)
            obj.ndofs = length(obj.vertices);
        end
        
        function computeShapeFunctions(obj)
            syms x y a b1 b2
            shapeFunc = cell(obj.ndofs,1);
            ndof = obj.ndofs;
            
            LHS = obj.applyLinearForm();
            for i = 1:ndof
                RHS = obj.computeLinearFormValues(i);
                c = obj.computeShapeFunctionCoefficients(LHS,RHS);
                shapeFunc{i} = matlabFunction(c.a+c.b1*x+c.b2*y,'Vars',[x y]);
            end
            obj.shapeFunctions = shapeFunc;
        end
        
        function c = computeShapeFunctionCoefficients(~,LHS,RHS)
            syms a b1 b2
            c = solve(LHS' == RHS,[a b1 b2]);
        end
        
        function RHS = computeLinearFormValues(obj,i)
            RHS = zeros(obj.ndofs,1);
            RHS(i) = 1;
        end
        
        function LHS = applyLinearForm(obj)
            syms x y a b1 b2
            baseShapeFunc = a+b1*x+b2*y;
            
            LHS(1) = obj.lineIntegral(baseShapeFunc,obj.vertices(2,:),obj.vertices(3,:));
            LHS(2) = obj.lineIntegral(baseShapeFunc,obj.vertices(1,:),obj.vertices(3,:));
            LHS(3) = obj.lineIntegral(baseShapeFunc,obj.vertices(2,:),obj.vertices(1,:));
        end
        
        function F = lineIntegral(~,func,pointA,pointB)
            syms x y t real
            x1 = pointA(1); y1 = pointA(2);
            x2 = pointB(1); y2 = pointB(2);
            func = subs(func,x,x1 + t*(x2-x1));
            func = subs(func,y,y1 + t*(y2-y1));
            F = int(func,t,0,1);
        end
        
    end
    
end