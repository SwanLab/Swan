classdef NedelecElement3D < handle
   
    properties (Access = private)
       domainK
        xSym
       ySym
       zSym
       a1Sym
       a2Sym
       a3Sym
       b1Sym
       b2Sym
       b3Sym
       shapeFunctionsSym
    end
    
    
    properties (Access = public)
        ndofs
        shapeFunctions
    end
    
    
    methods (Access = public)
    
        function obj = NedelecElement3D()
            obj.init();
        end
        
        function plotShapeFunctions(obj)
            nodes = [0,0,0;0,1/3,0;0,2/3,0;0,1,0;1/3,0,0;1/3,1/3,0;1/3,2/3,0;2/3,0,0;2/3,1/3,0;2/3,1/3,0;1,0,0;
                     0,0,1/3;0,1/3,1/3;0,2/3,1/3;1/3,0,1/3;1/3,1/3,1/3;2/3,0,1/3;0,0,2/3;0,1/3,2/3;1/3,0,2/3;0,0,1];
            m = obj.createPlotMesh();
            for s = 1:6
                figure();
                m.plot();
                x = zeros(length(nodes),3);
                for i = 1:length(nodes)
                    x(i,:) = obj.shapeFunctions{s}(nodes(i,1),nodes(i,2),nodes(i,3));
                end
                quiver3(nodes(:,1),nodes(:,2),nodes(:,3),x(:,1),x(:,2),x(:,3));
                xlabel('x'); ylabel('y'); zlabel('z')
                title("Shape Function (s = "+string(s-1)+")");
                grid on
                hold off
            end
        end
    
    end
    
    
    methods (Access = private)
       
        function init(obj)
            obj.domainK = Simplicial3D();
            
            obj.xSym = sym('x','real');
            obj.ySym = sym('y','real');
            obj.zSym = sym('z','real');
            obj.a1Sym = sym('a1','real');
            obj.a2Sym = sym('a2','real');
            obj.a3Sym = sym('a3','real');
            obj.b1Sym = sym('b1','real');
            obj.b2Sym = sym('b2','real');
            obj.b3Sym = sym('b3','real');
            
            obj.computeNdof();
            obj.computeShapeFunctionsSym();
            obj.computeShapeFunctions();
        end
        
        function computeNdof(obj)
            obj.ndofs = length(obj.domainK.tangentVectors);
        end
                
        function F = lineIntegral(obj,func,pointA,pointB)
            x = obj.xSym; y = obj.ySym; z = obj.zSym;
            t = sym('t','real');
            x1 = pointA(1); y1 = pointA(2); z1 = pointA(3);
            x2 = pointB(1); y2 = pointB(2); z2 = pointB(3);
            
            func = subs(func,x,x1 + t*(x2-x1));
            func = subs(func,y,y1 + t*(y2-y1));
            func = subs(func,z,z1 + t*(z2-z1));
            F = int(func,t,0,1);
        end
        
        function computeShapeFunctionsSym(obj)
            x = obj.xSym; y = obj.ySym; z = obj.zSym;
            ndof = obj.ndofs;
            shFunc = cell(ndof,1);
            
            LHS = obj.applyLinearForms();
            for s = 1:ndof
                RHS = obj.computeLinearFormValues(s);
                c = obj.computeShapeFunctionCoefficients(LHS,RHS);
                shFunc{s} = [c.a1+c.b3*y-c.b2*z,c.a2+c.b1*z-c.b3*x,c.a3+c.b2*x-c.b1*y];
            end
            
            obj.shapeFunctionsSym = shFunc;
        end
        
        function computeShapeFunctions(obj)
            ndof = obj.ndofs;
            shFunc = cell(ndof,1);
            shFuncSym = obj.shapeFunctionsSym;
            x = obj.xSym; y = obj.ySym; z = obj.zSym;
            
            for s = 1:ndof
                shFunc{s} = matlabFunction(shFuncSym{s},'Vars',[x y z]);
            end
            
            obj.shapeFunctions = shFunc;
        end
        
        function LHS = applyLinearForms(obj)
            x = obj.xSym; y = obj.ySym; z = obj.zSym;
            a1 = obj.a1Sym; a2 = obj.a2Sym; a3 = obj.a3Sym;
            b1 = obj.b1Sym; b2 = obj.b2Sym; b3 = obj.b3Sym;
            v = obj.domainK.vertices;
            ndof = obj.ndofs;
            t = obj.domainK.tangentVectors;
            
            N = [a1+b3*y-b2*z,a2+b1*z-b3*x,a3+b2*x-b1*y];
            for s = 1:ndof
                Nt(s) = dot(N,t(s,:));
            end
            
            LHS(1) = obj.lineIntegral(Nt(1),v(2,:),v(3,:));
            LHS(2) = obj.lineIntegral(Nt(2),v(3,:),v(1,:));
            LHS(3) = obj.lineIntegral(Nt(3),v(2,:),v(1,:));
            LHS(4) = obj.lineIntegral(Nt(4),v(1,:),v(4,:));
            LHS(5) = obj.lineIntegral(Nt(5),v(2,:),v(4,:));
            LHS(6) = obj.lineIntegral(Nt(6),v(3,:),v(4,:));
        end
        
        function RHS = computeLinearFormValues(obj,s)
            ndof = obj.ndofs;
            
            RHS = zeros(1,ndof);
            RHS(s) = obj.domainK.edgesLength(s);
        end
        
        function c = computeShapeFunctionCoefficients(obj,LHS,RHS)
            a1 = obj.a1Sym; a2 = obj.a2Sym; a3 = obj.a3Sym;
            b1 = obj.b1Sym; b2 = obj.b2Sym; b3 = obj.b3Sym;
            
            eq = LHS == RHS;
            c = solve(eq,[a1 a2 a3 b1 b2 b3]);
        end
        
        function m = createPlotMesh(obj)
            s.coord = obj.domainK.vertices;
            s.connec = obj.domainK.connectivities;
            
            m = Mesh(s);
        end
        
    end
    
end