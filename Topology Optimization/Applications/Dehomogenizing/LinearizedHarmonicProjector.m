classdef LinearizedHarmonicProjector < handle
    
    properties (Access = private)
        internalDOFs
        massMatrix
        stiffnessMatrix
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
            obj.createInternalDOFs();
            obj.computeMassMatrix();
            obj.computeStiffnessMatrix();
            obj.createSolver()
        end
        
        function lambda0 = computeInitalLambda(obj)
            iDOFs   = obj.internalDOFs;
            lambda0 = zeros(length(iDOFs),1);            
        end

        function eta0 = computeInitialEta(obj)
            npnod  = obj.mesh.nnodes;
            eta0    = zeros(npnod,1);
        end

        function [bNew,lambda] = solveProblem(obj,rho,bBar,b0)            
            b = obj.projectUnitBall(squeeze(b0.fValues)');
            i = 1;
            isErrorLarge = true;
            lambda    = obj.computeInitalLambda();
            eta       = obj.computeInitialEta();
            etaOld    = eta;
            lambdaOld = lambda;
            rhs      = obj.computeRHS(rho,bBar);
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
            iDOFs = obj.internalDOFs;
            resDT(iDOFs) = resD;
            figure(59)
            clf
            trisurf(obj.mesh.connec,x,y,resDT)
            view(0,90)
            colorbar            
        end

        function [iX,iY,iL,iE] = computeIndex(obj)
            f = P1Function.create(obj.mesh, 1);
            npnod = f.dim.ndofs;
            iDOFs = obj.internalDOFs;
            iX = 1:npnod;
            iY = (npnod) + (1:npnod);
            iL = (2*npnod + 1):(2*npnod +length(iDOFs));
            iE = (2*npnod + length(iDOFs) + 1):(2*npnod + length(iDOFs) + npnod);
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

        function createInternalDOFs(obj)
           b     = obj.boundaryMesh;
           iDOFs = setdiff(1:obj.mesh.nnodes,b); 
           obj.internalDOFs = iDOFs;
        end
        
        function computeMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = P1Function.create(obj.mesh, 1);
            s.trial = P1Function.create(obj.mesh, 1);
            s.quadratureOrder = 'QUADRATICMASS';
            lhs = LHSintegrator.create(s);
            M = lhs.compute();
            %M = diag(sum(M));
            %M = eye(size(M,1));
            obj.massMatrix = M;
        end

        function computeStiffnessMatrix(obj)        
            s.test  = P1Function.create(obj.mesh, 1);
            s.trial = P1Function.create(obj.mesh, 1);
            s.mesh         = obj.mesh;
            s.type         = 'StiffnessMatrix';
            lhs = LHSintegrator.create(s);
            K = lhs.compute();
            obj.stiffnessMatrix = K;
        end

        function [CX,CY,DX,DY,EX,EY,Kxx,Kxy,Kyy,Mxx,Mxy,Myy,Mrho] = computeAdvectionMatrix(obj,rho,vH)
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
%            s.field        = obj.field;
            s.type         = 'AdvectionMatrix';
            s.b            = vH;
            s.rho          = rho;
            s.quadType     = 'CUBIC';
            
            lhs = LHSintegrator.create(s);
            [CX,CY,DX,DY,EX,EY,Kxx,Kxy,Kyy,Mxx,Mxy,Myy,Mrho] = lhs.compute();
       end        

       function Ared = computeReducedAdvectionMatrix(obj,A)
           iDOFs = obj.internalDOFs;
           Ared = A(:,iDOFs);
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
           K    = obj.stiffnessMatrix;
           Zb   = zeros(size(M));
           Zbred = obj.computeReducedAdvectionMatrix(Zb);
           Z    = obj.computeZeroMatrix();

           %Ex = diag(sum(Ex));
           %Ey = diag(sum(Ey));
           Ex2 = diag(vH(:,1));
           Ey2 = diag(vH(:,2));
           l = 1;
           l2 = 0;
           lhs  = [Mrho+l*(Kxx+Mxx)+l2*K,Zb-l*(Kxy+Mxy),(-Dx+Cx),Ex;...
                   Zb-l*(Kxy+Mxy),Mrho+l*(Kyy+Myy)+l2*K,(Dy-Cy) ,Ey;...
                  (-Dx+Cx)',(Dy-Cy)',Z,Zbred';...
                   Ex2',Ey2',Zbred,Zb];
           %lhs  = [M,Zb,(-Dx+Cx),Ex;Zb,M,(Dy-Cy),Ey;(-Dx+Cx)',(Dy-Cy)',Z,Zbred';Ex',Ey',Zbred,Zb];
           obj.LHS = lhs;
       end

       function Z = computeZeroMatrix(obj)
           iDOFs = obj.internalDOFs;
           Z     = zeros(length(iDOFs),length(iDOFs));           
       end

        function rhsV = computeRHS(obj,rho,bBar)
            bBarV = squeeze(bBar.fValues)';

            % RHS not working for P0
%             s.fValues = bBarV(:,1);
%             s.mesh     = obj.mesh;
%             b1 = P0Function(s);
% 
%             s.fValues = bBarV(:,2);
%             s.mesh     = obj.mesh;
%             b2 = P0Function(s);
% 
%              s.mesh = obj.mesh;
%              s.fValues = zeros(obj.mesh.nnodes,1);
%              test = P1Function(s);
%              s.mesh  = obj.mesh;
%              s.type = 'ShapeFunction';
%              rhs    = RHSintegrator.create(s);
%              rhs1   = rhs.compute(b1,test);
%              rhs2   = rhs.compute(b2,test);


            bBar1 = bBar.project('P1');
            M = obj.massMatrix;
            w = 4*(1-rho).*rho;
            rhs1 = M*(bBar1.fValues(:,1).*w);
            rhs2 = M*(bBar1.fValues(:,2).*w);


            iDOFs = obj.internalDOFs;
            Z     = zeros(length(iDOFs),1);
            I     = ones(size(M,1),1);

%            rhs = [rhs1;rhs2;Z;M*I];
            %M2 = diag(sum(M));
           % rhs = [rhs1;rhs2;Z;M*I];
            rhsV = [rhs1;rhs2;Z;I];
            

        end

      

    end
    
end