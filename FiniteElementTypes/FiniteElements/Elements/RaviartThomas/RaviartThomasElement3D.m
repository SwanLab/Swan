classdef RaviartThomasElement3D < handle
   
    properties (Access = private)
       domainK
       shapeFunctionsSym
       xSym
       ySym
       zSym
       a1Sym
       a2Sym
       a3Sym
       bSym
    end
    
    
    properties (Access = public)
        shapeFunctions
        ndofs
    end
    
    
    methods (Access = public)
    
        function obj = RaviartThomasElement3D()
            obj.init();
        end
        
        function plotShapeFunctions(obj)
            set(groot,'defaulttextinterpreter','latex');
            ndof = obj.ndofs;
            nodes = [0,0,0;0,1/3,0;0,2/3,0;0,1,0;1/3,0,0;1/3,1/3,0;1/3,2/3,0;2/3,0,0;2/3,1/3,0;2/3,1/3,0;1,0,0;0,0,1/3;0,1/3,1/3;0,2/3,1/3;1/3,0,1/3;1/3,1/3,1/3;2/3,0,1/3;0,0,2/3;0,1/3,2/3;1/3,0,2/3;0,0,1];
            m = obj.createPlotMesh();
            for s = 1:ndof
                figure();
                m.plot();
                X = obj.shapeFunctions{s}(nodes(:,1),nodes(:,2),nodes(:,3));
                quiver3(nodes(:,1),nodes(:,2),nodes(:,3),double(X(:,1)),double(X(:,2)),double(X(:,3)),'k');
                xlabel('x'); ylabel('y'); zlabel('z');
                title("Shape Function (s = "+string(s-1)+")");
                grid on
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
            obj.bSym = sym('b','real');
            
            obj.computeNdof();
            obj.computeShapeFunctionsSym();
            obj.computeShapeFunctions();
        end
        
        function computeNdof(obj)
            obj.ndofs = length(obj.domainK.vertices);
        end
        
        function F = planeIntegral(obj,func,pointA,pointB,pointC)
            x = obj.xSym;
            y = obj.ySym;
            z = obj.zSym;
            u = sym('u','real');
            v = sym('v','real');
            
            paramx = (1-u-v)*pointA(1) + u*pointB(1) + v*pointC(1);
            paramy = (1-u-v)*pointA(2) + u*pointB(2) + v*pointC(2);
            paramz = (1-u-v)*pointA(3) + u*pointB(3) + v*pointC(3);
            func = subs(func,x,paramx);
            func = subs(func,y,paramy);
            func = subs(func,z,paramz);

            F = int(int(func,u,0,1),v,0,1);
        end
        
        function computeShapeFunctionsSym(obj)
            x = obj.xSym;
            y = obj.ySym;
            z = obj.zSym;
            ndof = obj.ndofs;
            shapeFunc = cell(ndof,1);
            
            LHS = obj.applyLinearForm();
            for s = 1:ndof
                RHS = obj.computeLinearFormValues(s);
                c = obj.computeShapeFunctionCoefficients(LHS,RHS);
                shapeFunc{s} = [c.a1+c.b*x,c.a2+c.b*y,c.a3+c.b*z];
            end
            
            obj.shapeFunctionsSym = shapeFunc;
        end
        
        function computeShapeFunctions(obj)
            ndof = obj.ndofs;
            shFunc = cell(ndof,1);
            shFuncSym = obj.shapeFunctionsSym;
            x = obj.xSym;
            y = obj.ySym;
            z = obj.zSym;
            
            for s = 1:ndof
                shFunc{s} = matlabFunction(shFuncSym{s},'Vars',[x y z]);
            end
            
            obj.shapeFunctions = shFunc;
        end
        
        function c = computeShapeFunctionCoefficients(obj,LHS,RHS)
            a1 = obj.a1Sym;
            a2 = obj.a2Sym;
            a3 = obj.a3Sym;
            b = obj.bSym;
            
            c = solve(LHS == RHS,[a1 a2 a3 b]);
        end
        
        function LHS = applyLinearForm(obj)
            x = obj.xSym;
            y = obj.ySym;
            z = obj.zSym;
            a1 = obj.a1Sym;
            a2 = obj.a2Sym;
            a3 = obj.a3Sym;
            b = obj.bSym;
            n = obj.domainK.normalVectors;
            v = obj.domainK.vertices;
            
            N = [a1+b*x,a2+b*y,a3+b*z];
            for j = 1:obj.ndofs
                Nn(j) = dot(N,n(j,:));
            end

            LHS(1) = obj.planeIntegral(Nn(1),v(2,:),v(3,:),v(4,:));
            LHS(2) = obj.planeIntegral(Nn(2),v(3,:),v(1,:),v(4,:));
            LHS(3) = obj.planeIntegral(Nn(3),v(1,:),v(2,:),v(4,:));
            LHS(4) = obj.planeIntegral(Nn(4),v(1,:),v(2,:),v(3,:));
        end
        
        function RHS = computeLinearFormValues(obj,s)
            ndof = obj.ndofs;
            RHS = zeros(1,ndof);
            
            RHS(s) = obj.domainK.facesArea(s);
        end
        
        function m = createPlotMesh(obj)
            s.coord = obj.domainK.vertices;
            s.connec = obj.domainK.connectivities;
            
            m = Mesh(s);
        end
        
    end
    
end