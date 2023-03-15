classdef NedelecElement2D < handle
   
    properties (Access = private)
       domainK
       xSym
       ySym
       a1Sym
       a2Sym
       bSym
       shapeFunctionsSym
    end
    
    
    properties (Access = public)
        ndofs
        shapeFunctions
    end
    
    
    methods (Access = public)
    
        function obj = NedelecElement2D()
            obj.init();
        end
        
        function plotShapeFunctions(obj)
            set(groot,'defaulttextinterpreter','latex');
            figure();
            m = obj.createPlotMeshR();
            mm = obj.createPlotMesh();
            for s = 1:3
                subplot(1,3,s)
                mm.plot();
                x = obj.shapeFunctions{s}(m.coord(:,1),m.coord(:,2));
                quiver(m.coord(:,1),m.coord(:,2),x(:,1),x(:,2),'k');
                title("Shape function (s = "+string(s-1)+")");
                xlabel('x'); ylabel('y');
                xlim([-0.5 1.5]); ylim([-0.5 1.5]);
                grid on
            end
        end
    
    end
    
    
    methods (Access = private)
       
        function init(obj)
            obj.domainK = Simplicial2D();
            obj.computeNdof();
            obj.computeShapeFunctionsSym();
            obj.computeShapeFunctions();
        end
        
        function computeNdof(obj)
            obj.ndofs = length(obj.domainK.vertices);
        end
        
        function computeShapeFunctionsSym(obj)
            obj.xSym = sym('x','real');
            obj.ySym = sym('y','real');
            x = obj.xSym;
            y = obj.ySym;
            ndof = obj.ndofs;
            shFunc = cell(ndof,1);
            
            LHS = obj.applyLinearForm();
            for s = 1:ndof
                RHS = obj.computeLinearFormValues(s);
                c = obj.computeShapeFunctionCoefficients(LHS,RHS);
                shFunc{s} = [c.a1-c.b*y,c.a2+c.b*x];
            end
            
            obj.shapeFunctionsSym = shFunc;
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
        
        function c = computeShapeFunctionCoefficients(obj,LHS,RHS)
            a1 = obj.a1Sym; a2 = obj.a2Sym; b = obj.bSym;
            eq = LHS == RHS;
            c = solve(eq,[a1 a2 b]);
        end
        
        function LHS = applyLinearForm(obj)
            obj.a1Sym = sym('a1','real');
            obj.a2Sym = sym('a2','real');
            obj.bSym = sym('b','real');
            x = obj.xSym; y = obj.ySym;
            a1 = obj.a1Sym; a2 = obj.a2Sym; b = obj.bSym;
            ndof = obj.ndofs;
            
            N = [a1-b*y,a2+b*x];
            for s = 1:ndof
                Nt(s) = dot(N,obj.domainK.tangentVectors(s,:));
            end
            
            v = obj.domainK.vertices;
            LHS(1) = obj.lineIntegral(Nt(1),v(2,:),v(3,:));
            LHS(2) = obj.lineIntegral(Nt(2),v(3,:),v(1,:));
            LHS(3) = obj.lineIntegral(Nt(3),v(1,:),v(2,:));
        end
        
        function RHS = computeLinearFormValues(obj,i)
            RHS = zeros(1,obj.ndofs);
            RHS(i) = 1;
        end
        
        function m = createPlotMeshR(obj)
            s.coord = obj.domainK.vertices;
            s.connec = [1 2 3];
            m = Mesh(s);
            for i=1:3
                m = m.remesh(2);
            end
        end
        
        function m = createPlotMesh(obj)
            s.coord = obj.domainK.vertices;
            s.connec = [1 2 3];
            m = Mesh(s);
        end
        
    end
    
end