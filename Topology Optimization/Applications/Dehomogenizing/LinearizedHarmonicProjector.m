classdef LinearizedHarmonicProjector < handle
    
    properties (Access = private)
        dim
        massMatrix
        advectionMatrixX
        advectionMatrixY
        reducedAdvectionMatrixX
        reducedAdvectionMatrixY
        LHS
        solver
    end    

    properties (Access = private)
       mesh   
       boundaryMesh
    end
    
    methods (Access = public)
        
        function obj = LinearizedHarmonicProjector(cParams)
            obj.init(cParams);            
            obj.createDimension();
            obj.computeMassMatrix();
            obj.createSolver()
        end
        
        function lambda0 = computeInitalLambda(obj)
            b    = obj.boundaryMesh;
            nInt = setdiff(1:obj.dim.npnod,b);
            lambda0 = zeros(length(nInt),1);
        end

        function [v,lambda] = solveProblem(obj,v0,vH)            
            vH = obj.projectUnitBall(vH);
            obj.plotOrientation(vH,1);
            i = 1;
            isErrorLarge = true;
            lambda = obj.computeInitalLambda();
            tau =1000;
            while isErrorLarge
                obj.plotOrientation(vH,2);

                obj.computeLHS(vH);
                lhs    = obj.LHS;
                rhs    = obj.computeRHS(v0);


                % gradient                            
                I = 1:2*obj.dim.npnod;
                Il = (2*obj.dim.npnod + 1):length(rhs);
                x = [vH(:,1);vH(:,2);lambda];
                res  = lhs*x - rhs;
                resP = res(I);
                resD = res(Il);


                b    = obj.boundaryMesh;
                nInt = setdiff(1:obj.dim.npnod,b);                               
                x = obj.mesh.coord(:,1);
                y = obj.mesh.coord(:,2);
                z = zeros(size(x));
                z(nInt) = resD;
                figure()
                trisurf(obj.mesh.connec,x,y,z)
                view(0,90)
                colorbar

%                 sol = x - tau*(res);
%                 lambda = lambda + (lhs(Il,I)*x );
%                 indexX = 1:obj.dim.npnod;
%                 indexY = (obj.dim.npnod) + (1:obj.dim.npnod);                
%                 v(:,1) = sol(indexX,1);
%                 v(:,2) = sol(indexY,1);                

                % Picard
                sol    = obj.solver.solve(lhs,rhs); 
                indexX = 1:obj.dim.npnod;
                indexY = (obj.dim.npnod) + (1:obj.dim.npnod);
                indexL = (2*(obj.dim.npnod) + 1):length(sol);
                v(:,1) = sol(indexX,1);
                v(:,2) = sol(indexY,1);
                lambda = sol(indexL,1);

                
                
                % err(i) = norm(vH(:)-v(:))/norm(v(:));
                err(i)  = norm(res);
                errP(i) = norm(resP);                
                errD(i) = norm(resD);

                isErrorLarge = err(i) > 1e-12;
                i = i + 1;                
                figure(99)
                plot(log10(err)) 
                theta = 0.99;
                v = obj.projectUnitBall(v);
                vH = theta*v + (1-theta)*vH ;     
                vH = obj.projectUnitBall(vH);
            end



        end

        function tp = projectUnitBall(obj,t)
            uP = UnitBallProjector([]);
            tp = uP.project(t);
        end

        function plotOrientation(obj,t,iFigure)
            figure(100)
            subplot(2,2,iFigure)            
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            tx = t(:,1);
            ty = t(:,2);
            quiver(x,y,tx,ty);
            tp = obj.projectUnitBall(t);
            tx = tp(:,1);
            ty = tp(:,2);
            subplot(2,2,2*(2-1)+iFigure)            
            quiver(x,y,tx,ty);            
        end

        function optPrim = computePrimalOptimaility(obj,lambda,v,vH)
            M        = obj.massMatrix;  
            Ared     = obj.reducedAdvectionMatrix;  
            gradPrim = M*(v-vH) + (Ared'*lambda);
            optPrim  = norm(gradPrim);
        end

        function vH = projectByDual(obj,v)
            Kred = obj.reducedAdvectionMatrix;            
            M    = obj.massMatrix;            
            lhs = Kred*(M\Kred');
            rhs = Kred*v;            
            lambda  = obj.solver.solve(lhs,rhs);
            vH = v - M\(Kred'*lambda);
        end    

        function cost = computeCost(obj,v,vH)
            M    = obj.massMatrix;     
            cost = 0;
            for i = 1:size(v,2)
                vI   = v(:,i);
                incV = vI- vH(:,i);
                cost = cost + incV'*M*incV/(vI'*M*vI);            
            end
        end

        function optDual = computeDualOptimality(obj,vH)
            Ared = obj.reducedAdvectionMatrix;            
            optDual = norm(Ared*vH);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.mesh             = cParams.mesh;
           obj.boundaryMesh     = cParams.boundaryMesh;
        end
        
        function createDimension(obj)
            q = Quadrature();
            q = q.set(obj.mesh.type);
            s.mesh = obj.mesh;
            s.pdim = '1D';
            s.ngaus = q.ngaus;
            d = DimensionVariables(s);
            d.compute();
            obj.dim = d;
        end        
        
        function computeMassMatrix(obj)
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.npnod        = obj.mesh.npnod;
            s.type         = 'MassMatrix';
            s.dim          = obj.dim;
            s.quadType     = 'QUADRATIC';
         
            lhs = LHSintegrator.create(s);
            M = lhs.compute();
            %M = diag(sum(M));
            %M = eye(size(M,1));
            obj.massMatrix = M;
        end

        function [CX,CY,DX,DY] = computeAdvectionMatrix(obj,vH)
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.npnod        = obj.mesh.npnod;
            s.type         = 'AdvectionMatrix';
            s.dim          = obj.dim;
            s.b            = vH;
            lhs = LHSintegrator.create(s);
            [CX,CY,DX,DY] = lhs.compute();
       end        

       function Ared = computeReducedAdvectionMatrix(obj,A)
           b    = obj.boundaryMesh;
           nInt = setdiff(1:obj.dim.npnod,b);
           Ared = A(:,nInt);
       end

       function createSolver(obj)
            s = Solver.create();
            obj.solver = s;
       end


       function  computeLHS(obj,vH)
           [Cx,Cy,Dx,Dy] = obj.computeAdvectionMatrix(vH);
           Cx = obj.computeReducedAdvectionMatrix(Cx);
           Cy = obj.computeReducedAdvectionMatrix(Cy);
           Dx = obj.computeReducedAdvectionMatrix(Dx);
           Dy = obj.computeReducedAdvectionMatrix(Dy);
           M    = obj.massMatrix;
           Zb   = zeros(size(M));
           Z    = obj.computeZeroFunction();
           lhs  = [M,Zb,(0*Cx+1*Dx);Zb,M,(0*Cy+1*Dy);(0*Cx + Dx)',(0*Cy + Dy)',Z];
           obj.LHS = lhs;
       end

       function Z = computeZeroFunction(obj)
           b    = obj.boundaryMesh;
           nInt = setdiff(1:obj.dim.npnod,b);           
           Z    = zeros(length(nInt),length(nInt));           
       end

        function rhs = computeRHS(obj,v)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            s.mesh  = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.npnod = obj.dim.npnod;
            s.dim   = obj.dim;
            s.type = 'SIMPLE';
            int = Integrator.create(s);
            rhs1 = int.integrateFnodal(v(:,1),q.order);
            rhs2 = int.integrateFnodal(v(:,2),q.order);
            b = obj.boundaryMesh;
            nInt = setdiff(1:obj.dim.npnod,b);
            Z   = zeros(length(nInt),1);


           [Cx,Cy,Dx,Dy] = obj.computeAdvectionMatrix(v);
           Cx = obj.computeReducedAdvectionMatrix(Cx);
           Cy = obj.computeReducedAdvectionMatrix(Cy);            
           res = 0*Cx'*v(:,1) + 0*Cy'*v(:,2);
            rhs = [rhs1;rhs2;Z+res];
        end

      

    end
    
end