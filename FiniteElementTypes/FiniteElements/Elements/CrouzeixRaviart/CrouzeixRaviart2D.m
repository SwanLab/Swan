classdef CrouzeixRaviart2D < handle
   
    properties (Access = private)
       xSym
       ySym
       aSym
       b1Sym
       b2Sym
       shapeFunctionsSym
       domainK
    end
    
    
    properties (Access = public)
        shapeFunctions
        ndofs
    end
    
    
    methods (Access = public)
    
        function obj = CrouzeixRaviart2D()
            obj.init();
        end
        
        function plotShapeFunctions(obj)
            set(groot,'defaulttextinterpreter','latex');
            ndof = obj.ndofs;
            
            m = obj.createPlotMesh();
            for s = 1:ndof
                figure();
                m.plot();
                trisurf(m.connec,m.coord(:,1),m.coord(:,2),obj.shapeFunctions{s}(m.coord(:,1),m.coord(:,2)));
                xlim([0 1]); ylim([0 1]);
                xlabel('x'); ylabel('y'); zlabel('z');
                title("Shape Function (s = "+string(s-1)+")");
                grid on
            end
        end
    
    end
    
    
    methods (Access = private)
       
        function init(obj)
            obj.domainK = Simplicial2D();
            
            obj.xSym = sym('x','real');
            obj.ySym = sym('y','real');
            obj.aSym = sym('a','real');
            obj.b1Sym = sym('b1','real');
            obj.b2Sym = sym('b2','real');
            
            obj.computeNdof();
            obj.computeShapeFunctionsSym();
            obj.computeShapeFunctions();
        end
        
        function computeNdof(obj)
            obj.ndofs = length(obj.domainK.vertices);
        end
        
        function computeShapeFunctionsSym(obj)
            x = obj.xSym;
            y = obj.ySym;
            shapeFunc = cell(obj.ndofs,1);
            ndof = obj.ndofs;
            
            LHS = obj.applyLinearForm();
            for s = 1:ndof
                RHS = obj.computeLinearFormValues(s);
                c = obj.computeShapeFunctionCoefficients(LHS,RHS);
                shapeFunc{s} = c.a+c.b1*x+c.b2*y;
            end
            
            obj.shapeFunctionsSym = shapeFunc;
        end
        
        function computeShapeFunctions(obj)
            ndof = obj.ndofs;
            shFunc = cell(ndof,1);
            shFuncSym = obj.shapeFunctionsSym;
            x = obj.xSym;
            y = obj.ySym;
            
            for s = 1:ndof
                shFunc{s} = matlabFunction(shFuncSym{s},'Vars',[x y]);
            end
            
            obj.shapeFunctions = shFunc;
        end
        
        function c = computeShapeFunctionCoefficients(obj,LHS,RHS)
            a = obj.aSym;
            b1 = obj.b1Sym;
            b2 = obj.b2Sym;
            
            c = solve(LHS' == RHS,[a b1 b2]);
        end
        
        function RHS = computeLinearFormValues(obj,s)
            RHS = zeros(obj.ndofs,1);
            RHS(s) = 1;
        end
        
        function LHS = applyLinearForm(obj)
            x = obj.xSym;
            y = obj.ySym;
            a = obj.aSym;
            b1 = obj.b1Sym;
            b2 = obj.b2Sym;
            N = a+b1*x+b2*y;
            v = obj.domainK.vertices;
            
            LHS(1) = obj.lineIntegral(N,v(2,:),v(3,:));
            LHS(2) = obj.lineIntegral(N,v(1,:),v(3,:));
            LHS(3) = obj.lineIntegral(N,v(2,:),v(1,:));
        end
        
        function F = lineIntegral(obj,func,pointA,pointB)
            x = obj.xSym;
            y = obj.ySym;
            t = sym('t','real');
            x1 = pointA(1); y1 = pointA(2);
            x2 = pointB(1); y2 = pointB(2);
            
            func = subs(func,x,x1 + t*(x2-x1));
            func = subs(func,y,y1 + t*(y2-y1));
            F = int(func,t,0,1);
        end
        
        function m = createPlotMesh(obj)
            s.coord = obj.domainK.vertices;
            s.connec = [1 2 3];
            m = Mesh(s);
        end
        
    end
    
end