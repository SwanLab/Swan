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
            lambda0 = 0*ones(length(nInt),1);
        end

        function [v,lambda] = solveProblem(obj,v0,vH)            
            vH = obj.projectUnitBall(vH);
            obj.plotOrientation(vH,1);
            i = 1;
            isErrorLarge = true;
            lambda = obj.computeInitalLambda();
            npnod = obj.dim.npnod;
            eta    = zeros(npnod,1);
            etaOld = eta;
            lambdaOld = lambda;
            tau = 1000;
            rhs    = obj.computeRHS(v0);
            t = 1;
            while isErrorLarge
                 if mod(i,10) == 0
                     obj.plotOrientation(vH,2);
                 end


                obj.computeLHS(vH);
                lhs    = obj.LHS;


                % gradient                            
                I = 1:2*npnod;
                nInt = obj.computeNint();
                Il =  (2*npnod + 1):(2*npnod +length(nInt));
                Ieta = (2*npnod + length(nInt) + 1):(2*npnod + length(nInt) + npnod);
                x = [vH(:,1);vH(:,2);lambda;eta];
                res  = lhs*x - rhs;
                resP = res(I);
                resD = res(Il);
                resE = res(Ieta);

%                 b    = obj.boundaryMesh;
%                 nInt = setdiff(1:obj.dim.npnod,b);                               
%                 x = obj.mesh.coord(:,1);
%                 y = obj.mesh.coord(:,2);
%                 z = zeros(size(x));
%                 z(nInt) = resD;
%                 figure()
%                 trisurf(obj.mesh.connec,x,y,z)
%                 view(0,90)
%                 colorbar

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
                indexL = (2*npnod + 1):(2*npnod +length(nInt));
                indexEta = (2*npnod + length(nInt) + 1):(2*npnod + length(nInt) + npnod);
                v(:,1) = sol(indexX,1);
                v(:,2) = sol(indexY,1);
                lambda = sol(indexL,1);
                eta    = sol(indexEta,1);

                
                
                % err(i) = norm(vH(:)-v(:))/norm(v(:));
                err(i)  = norm(res);
                errP(i) = norm(resP);                
                errD(i) = norm(resD);
                errE(i) = norm(resE);
                incX(i) = norm(v(:)-vH(:));
                incL(i) = norm(lambda-lambdaOld);
                incE(i) = norm(eta-etaOld);

                lambdaOld = lambda;
                etaOld = eta;
                isErrorLarge = err(i) > 1e-13;

                if mod(i,10) == 0
                figure(99)
                semilogy(1:i,([err;errP;errD;errE;incX;incL;incE])) 
                legend('res','resP','resD','resE','incX','incL','incE')
                end                
                                
                
                theta = 0.5;
             %   v = obj.projectUnitBall(v);
                vH = theta*v + (1-theta)*vH ;    
%                 if mod(t,1000) == 0
%                 vH = obj.projectUnitBall(vH);

              %  theta2 = 0;
              %  normvH = obj.computeNorm(vH);
               % normT = theta2*normvH + (1-theta2)*1;
              %  vH = vH./normT;



%                 t = t+1;
%                 end
                i = i + 1;
            end



        end

        function norm = computeNorm(obj,v)
         vx = v(:,1);
         vy = v(:,2);
         norm = sqrt(vx.^2 + vy.^2);
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

        function [CX,CY,DX,DY,EX,EY] = computeAdvectionMatrix(obj,vH)
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.npnod        = obj.mesh.npnod;
            s.type         = 'AdvectionMatrix';
            s.dim          = obj.dim;
            s.b            = vH;
            s.quadType     = 'QUADRATIC';
            
            lhs = LHSintegrator.create(s);
            [CX,CY,DX,DY,EX,EY] = lhs.compute();
       end        

       function Ared = computeReducedAdvectionMatrix(obj,A)
           nInt = obj.computeNint();
           Ared = A(:,nInt);
       end

       function createSolver(obj)
            s = Solver.create();
            obj.solver = s;
       end


       function  computeLHS(obj,vH)
           [Cx,Cy,Dx,Dy,Ex,Ey] = obj.computeAdvectionMatrix(vH);
           Cx = obj.computeReducedAdvectionMatrix(Cx);
           Cy = obj.computeReducedAdvectionMatrix(Cy);
           Dx = obj.computeReducedAdvectionMatrix(Dx);
           Dy = obj.computeReducedAdvectionMatrix(Dy);
           M    = obj.massMatrix;
           Zb   = zeros(size(M));
           Zbred = obj.computeReducedAdvectionMatrix(Zb);
           Z    = obj.computeZeroMatrix();
           lhs  = [M,Zb,(-Dx+Cx),Ex;Zb,M,(Dy-Cy),Ey;(-Dx+Cx)',(Dy-Cy)',Z,Zbred';Ex',Ey',Zbred,Zb];
           obj.LHS = lhs;
       end

       function Z = computeZeroMatrix(obj)
           nInt = obj.computeNint();
           Z    = zeros(length(nInt),length(nInt));           
       end

       function nInt = computeNint(obj)
            b    = obj.boundaryMesh;
           nInt = setdiff(1:obj.dim.npnod,b);                     
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
      %      rhs1 = int.integrateFnodal(v(:,1),q.order);
       %     rhs2 = int.integrateFnodal(v(:,2),q.order);


            M = obj.massMatrix;
            rhs1 = M*v(:,1);
            rhs2 = M*v(:,2);


            nInt = obj.computeNint();
            Z   = zeros(length(nInt),1);
            I   = ones(size(M,1),1);

            rhs = [rhs1;rhs2;Z;I];
        end

      

    end
    
end