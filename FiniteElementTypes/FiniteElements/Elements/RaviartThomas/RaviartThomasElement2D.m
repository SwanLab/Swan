classdef RaviartThomasElement2D < handle
   
    properties (Access = private)
       domainK
       shapeFunctionsSym
       shapeFunctionsDiffSym
       xSym
       ySym
       a1Sym
       a2Sym
       bSym
    end
    
    
    properties (Access = public)
        ndofs
        shapeFunctions
        shapeFunctionsDiff
        polynomialOrder
        type
    end
    
    
    methods (Access = public)
    
        function obj = RaviartThomasElement2D()
            obj.init();
            obj.polynomialOrder = 1; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
        
        function [z1,z2] = evaluate(obj,dofs,X,Y)
            coor = obj.computeNodesElement([X,Y]);
            z1 = zeros(length(X),1);
            z2 = zeros(length(X),1);
            
            for i = 1:obj.ndofs
                r = dofs(i).*obj.shapeFunctions{i}(coor(:,1),coor(:,2));
                z1 = z1 + r(:,1);
                z2 = z2 + r(:,2);
            end
            
            [z1,z2] = obj.computeVectMesh([X,Y],coor,[z1,z2]);
            
        end
        
        function [z1,z2] = computeVectMesh(obj,coords,coor,vects)
            vecpoint = coor+vects;
            base = coor(1:3,:);          
            M = (coords(1:3,:)-coords(1,:))'/base';
            N = coords(1,:)';
            v = (M*vecpoint'+N)-coords';
            z1 = v(1,:)';
            z2 = v(2,:)';
        end
        
        function coor = computeNodesElement(obj,coords)
            base = obj.domainK.vertices;          
            c = base(1:3,:);
            M = c'/(coords(1:3,:)-coords(1,:))';
            N = coords(1,:)';
            coor = (M*(coords-N')')';
        end
    
    end
    
    
    methods (Access = private)
       
        function init(obj)
            obj.domainK = Simplicial2D();
            
            obj.type = "RT";
            
            obj.xSym = sym('x','real');
            obj.ySym = sym('y','real');
            obj.a1Sym = sym('a1','real');
            obj.a2Sym = sym('a2','real');
            obj.bSym = sym('b','real');
            
            obj.computeNdof();
            obj.computeShapeFunctionsSym();
            obj.computeShapeFunctions();
            obj.computeShapeFunctionsDiffSym();
            obj.computeShapeFunctionsDiff();
        end
        
        function computeNdof(obj)
            obj.ndofs = length(obj.domainK.vertices);
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
        
        function computeShapeFunctionsSym(obj)
            x = obj.xSym;
            y = obj.ySym;
            ndof = obj.ndofs;
            shFunc = cell(ndof,1);
            
            LHS = obj.applyLinearForms();
            for s = 1:ndof
                RHS = obj.computeLinearFormsValues(s);
                c = obj.computeShapeFunctionCoefficients(LHS,RHS);
                shFunc{s} = [c.a1+c.b*x,c.a2+c.b*y];
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
        
        function c = computeShapeFunctionCoefficients(obj,LHS,RHS)
            a1 = obj.a1Sym;
            a2 = obj.a2Sym;
            b = obj.bSym;
            c = solve(LHS == RHS,[a1 a2 b]);
        end
        
        function LHS = applyLinearForms(obj,Nn)
            x = obj.xSym;
            y = obj.ySym;
            a1 = obj.a1Sym;
            a2 = obj.a2Sym;
            b = obj.bSym;
            v = obj.domainK.vertices;
            n = obj.domainK.normalVectors;
            
            N = [a1+b*x,a2+b*y];
            for s = 1:obj.ndofs
                Nn(s) = dot(N,n(s,:));
            end
            
            LHS(2) = obj.lineIntegral(Nn(3),v(2,:),v(3,:))*sqrt(2);
            LHS(3) = obj.lineIntegral(Nn(2),v(3,:),v(1,:));
            LHS(1) = obj.lineIntegral(Nn(1),v(1,:),v(2,:));
        end
        
        function RHS = computeLinearFormsValues(obj,s)
            ndof = obj.ndofs;
            
            RHS = zeros(1,ndof);
            RHS(s) = 1;
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
        
        function computeShapeFunctionsDiffSym(obj)
            f = obj.shapeFunctionsSym;
            shD = cell(length(f),2);
            ndof = obj.ndofs;
            x = obj.xSym;
            y = obj.ySym;
            
            for s = 1:ndof
                shD{s,1} = diff(f{s},x);
                shD{s,2} = diff(f{s},y);
            end
            
            obj.shapeFunctionsDiffSym = shD;
        end
        
        function computeShapeFunctionsDiff(obj)
            ndof = obj.ndofs;
            shD = cell(ndof,2);
            shDSym = obj.shapeFunctionsDiffSym;
            x = obj.xSym;
            y = obj.ySym;
            
            for s = 1:ndof
                shD{s,1} = matlabFunction(shDSym{s,1},'Vars',[x y]);
                shD{s,2} = matlabFunction(shDSym{s,2},'Vars',[x y]);
            end
            
            obj.shapeFunctionsDiff = shD;
        end
        
    end
    
end

