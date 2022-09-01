classdef LinearizedHarmonicProjector < handle
    
    properties (Access = private)
        field
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
            obj.createField();
            obj.computeMassMatrix();
            obj.createSolver()
        end
        
        function lambda0 = computeInitalLambda(obj)
            b    = obj.boundaryMesh;
            nInt = setdiff(1:obj.field.dim.ndofs,b);
            lambda0 = 0*ones(length(nInt),1);
        end

        function [bNew,lambda] = solveProblem(obj,rho,bBar,b0)            
            b = obj.projectUnitBall(b0);
            obj.plotOrientation(b,1);
            i = 1;
            isErrorLarge = true;
            lambda = obj.computeInitalLambda();
            npnod = obj.field.dim.ndofs;
            eta    = zeros(npnod,1);
            etaOld = eta;
            lambdaOld = lambda;
            rhs    = obj.computeRHS(rho,bBar);
            while isErrorLarge
             

                obj.computeLHS(rho,b);
                lhs    = obj.LHS;


                % gradient                            
                [iX,iY,iL,iE] = obj.computeIndex();
                x = [b(:,1);b(:,2);lambda;eta];
                res  = lhs*x - rhs;
                resP = res([iX,iY]);
                resD = res(iL);
                resE = res(iE);
           

                % Picard
                sol    = obj.solver.solve(lhs,rhs); 
                [iX,iY,iL,iE] = obj.computeIndex();
                bNew(:,1) = sol(iX,1);
                bNew(:,2) = sol(iY,1);
                lambda = sol(iL,1);
                eta    = sol(iE,1);

                
                
                err(i)  = norm(res);
                errP(i) = norm(resP);                
                errD(i) = norm(resD);
                errE(i) = norm(resE);
                incX(i) = norm(bNew(:)-b(:));
                incL(i) = norm(lambda-lambdaOld);
                incE(i) = norm(eta-etaOld);

                lambdaOld = lambda;
                etaOld = eta;
                isErrorLarge = err(i) > 1e-2;%1e-13;

                if mod(i,10) == 0
                    obj.plotOrientation(b,2);                    
                    figure(99)
                    semilogy(1:i,([err;errP;errD;errE;incX;incL;incE]))
                    legend('res','resP','resD','resE','incX','incL','incE')
                    errN = obj.computeErrorNorm(b);
                    obj.plotResiudalUnitBall(resE,errN)
                end                
                                
                
                theta = 0.11;%0.5;
             %   v = obj.projectUnitBall(v);
                b = theta*bNew + (1-theta)*b ;    


                i = i + 1;
            end



        end

        function plotHarmonicity(obj,resD)
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            resDT = zeros(size(x));
            nInt = obj.computeNint();
            resDT(nInt) = resD;
            figure(59)
            clf
            trisurf(obj.mesh.connec,x,y,resDT)
            view(0,90)
            colorbar            
        end

        function [iX,iY,iL,iE] = computeIndex(obj)
            npnod = obj.field.dim.ndofs;
            nInt  = obj.computeNint();
            iX = 1:npnod;
            iY = (npnod) + (1:npnod);
            iL = (2*npnod + 1):(2*npnod +length(nInt));
            iE = (2*npnod + length(nInt) + 1):(2*npnod + length(nInt) + npnod);
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

        function plotResiudalUnitBall(obj,resE,errN)
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            figure(55)
            clf
            subplot(1,2,1)
            trisurf(obj.mesh.connec,x,y,resE)
            view(0,90)
            colorbar
            subplot(1,2,2)
            trisurf(obj.mesh.connec,x,y,errN)     
            view(0,90)  
            colorbar
        end

        function plotOrientation(obj,t,iFigure)
            figure(100)
            subplot(2,2,iFigure)            
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            tx = t(:,1);
            ty = t(:,2);
           
            tp = obj.projectUnitBall(t);
            q = quiver(x,y,tx-tp(:,1),ty-tp(:,2),1.1);
            q.ShowArrowHead = 'off';
            tx = tp(:,1);
            ty = tp(:,2);
            subplot(2,2,2*(2-1)+iFigure)            
            quiver(x,y,tx,ty);            
        end

        function errN = computeErrorNorm(obj,t)
            tNorm = obj.computeNorm(t);
            errN = abs(1-tNorm);
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
             
        function createField(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1;
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'QUADRATIC';
            f = Field(s);
            obj.field = f;
        end
        
        function computeMassMatrix(obj)


            s.field        = obj.field;
            s.mesh         = obj.mesh;
            s.type         = 'MassMatrix';
            


         
            lhs = LHSintegrator.create(s);
            M = lhs.compute();
            %M = diag(sum(M));
            %M = eye(size(M,1));
            obj.massMatrix = M;
        end

        function [CX,CY,DX,DY,EX,EY,Kxx,Kxy,Kyy,Mxx,Mxy,Myy,Mrho] = computeAdvectionMatrix(obj,rho,vH)
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.field        = obj.field;
            s.type         = 'AdvectionMatrix';
            s.b            = vH;
            s.rho          = rho;
            s.quadType     = 'CUBIC';
            
            lhs = LHSintegrator.create(s);
            [CX,CY,DX,DY,EX,EY,Kxx,Kxy,Kyy,Mxx,Mxy,Myy,Mrho] = lhs.compute();
       end        

       function Ared = computeReducedAdvectionMatrix(obj,A)
           nInt = obj.computeNint();
           Ared = A(:,nInt);
       end

       function createSolver(obj)
            sS.type = 'DIRECT';
            s = Solver.create(sS);
            obj.solver = s;
       end


       function  computeLHS(obj,rho,vH)
          
           [Cx,Cy,Dx,Dy,Ex,Ey,Kxx,Kxy,Kyy,Mxx,Mxy,Myy,Mrho] = obj.computeAdvectionMatrix(rho,vH);
           Cx = obj.computeReducedAdvectionMatrix(Cx);
           Cy = obj.computeReducedAdvectionMatrix(Cy);
           Dx = obj.computeReducedAdvectionMatrix(Dx);
           Dy = obj.computeReducedAdvectionMatrix(Dy);
           M    = obj.massMatrix;
           Zb   = zeros(size(M));
           Zbred = obj.computeReducedAdvectionMatrix(Zb);
           Z    = obj.computeZeroMatrix();

           %Ex = diag(sum(Ex));
           %Ey = diag(sum(Ey));
           Ex2 = diag(vH(:,1));
           Ey2 = diag(vH(:,2));
           l = 1;
           lhs  = [Mrho+l*(Kxx+Mxx),Zb-l*(Kxy+Mxy)   ,(-Dx+Cx),Ex;...
                   Zb-l*(Kxy+Mxy),Mrho+l*(Kyy+Myy)   ,(Dy-Cy) ,Ey;...
                  (-Dx+Cx)',(Dy-Cy)',Z,Zbred';...
                   Ex2',Ey2',Zbred,Zb];
           %lhs  = [M,Zb,(-Dx+Cx),Ex;Zb,M,(Dy-Cy),Ey;(-Dx+Cx)',(Dy-Cy)',Z,Zbred';Ex',Ey',Zbred,Zb];
           obj.LHS = lhs;
       end

       function Z = computeZeroMatrix(obj)
           nInt = obj.computeNint();
           Z    = zeros(length(nInt),length(nInt));           
       end

       function nInt = computeNint(obj)
            b    = obj.boundaryMesh;
           nInt = setdiff(1:obj.field.dim.ndofs,b);                     
       end

        function rhs = computeRHS(obj,rho,v)
%             q = Quadrature.set(obj.mesh.type);
%             q.computeQuadrature('LINEAR');
%             s.mesh  = obj.mesh;
%             s.globalConnec = obj.mesh.connec;
%             s.npnod = obj.dim.ndofs;
%             s.dim   = obj.dim;
%             s.type = 'SIMPLE';
      %      int = Integrator.create(s);
      %      rhs1 = int.integrateFnodal(v(:,1),q.order);
       %     rhs2 = int.integrateFnodal(v(:,2),q.order);


            M = obj.massMatrix;
            w = 4*(1-rho).*rho;
            rhs1 = M*(v(:,1).*w);
            rhs2 = M*(v(:,2).*w);


            nInt = obj.computeNint();
            Z   = zeros(length(nInt),1);
            I   = ones(size(M,1),1);

%            rhs = [rhs1;rhs2;Z;M*I];
            %M2 = diag(sum(M));
           % rhs = [rhs1;rhs2;Z;M*I];
            rhs = [rhs1;rhs2;Z;I];
            

        end

      

    end
    
end