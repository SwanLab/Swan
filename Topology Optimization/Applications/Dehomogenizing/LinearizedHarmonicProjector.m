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
       l1
       l2
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

        function hRes = evaluateHarmonicResidual(obj,res,b)
            quad = Quadrature.set(obj.mesh.type);            
            quad.computeQuadrature('QUADRATICMASS');
            xV = quad.posgp;

            bG  = b.evaluate(xV);
            b1  = bG(1,:,:);
            b2  = bG(2,:,:);
            gbG = b.evaluateGradient(xV);
            gb1 = gbG.fValues(1:2,:,:);
            gb2 = gbG.fValues(3:4,:,:);

            s.fValues(1,:,:) = -b2.*gb1(1,:,:) + b1.*gb2(1,:,:);
            s.fValues(2,:,:) = -b2.*gb1(2,:,:) + b1.*gb2(2,:,:);
            s.mesh = obj.mesh;
            s.quadrature = xV;
            hRes = FGaussDiscontinuousFunction(s);
            
%             b2(1,:) = squeeze(b0.fValues(2,1,:));
%             b2(2,:) = squeeze(b0.fValues(2,1,:));


% 
%             s.fValues = atan2(b.fValues(:,2),b.fValues(:,1)); 
%             s.mesh = obj.mesh; 
%             beta = P1Function(s);
%             gradB = beta.computeGradient(quad);            
%             s.fValues = (gradB.fValues').^2;
%             s.mesh = obj.mesh; 
%             gH = P0Function(s);
%             gH.plot();
% 
%             gB = b.computeGradient(quad);
%             gBV = gB.fValues;
% 
%             s.fValues = gBV([1,2],:)';
%             s.mesh    = obj.mesh;
%             gb1 = P0Function(s);
% 
%             s.fValues = gBV([3,4],:)';
%             s.mesh    = obj.mesh;
%             gb2 = P0Function(s);
% 
%             b0 = b.project('P0');
%             b1(1,:) = squeeze(b0.fValues(1,1,:));
%             b1(2,:) = squeeze(b0.fValues(1,1,:));
%             b2(1,:) = squeeze(b0.fValues(2,1,:));
%             b2(2,:) = squeeze(b0.fValues(2,1,:));
% 
%             gradB2.fValues = -b1.*squeeze(gb2.fValues) + b2.*squeeze(gb1.fValues);
%             s.fValues = (gradB2.fValues').^2;
%             s.mesh = obj.mesh; 
%             gH2 = P0Function(s);
%             gH2.plot();       
% 
%             hRes = gH2;
% 
%             s.fValues = abs(squeeze(gradB.fValues)'-squeeze(gradB2.fValues)');s.mesh = obj.mesh; e = P0Function(s); e.plot            
% 
%             [iX,~,iL,~] = obj.computeIndex();
%             resH  = res(iL);
%             resHv = zeros(length(iX),1);
%             iDOFs = obj.internalDOFs;
%             resHv(iDOFs) = abs(resH);
%             s.fValues = resHv;
%             s.mesh    = obj.mesh;
%             hRes = P1Function(s);
        end

        function lRes = evaluateLossResidual(obj,res,bBar,b)
            bNorm1  = obj.L2norm(bBar.fValues(:,1));
            bNorm2  = obj.L2norm(bBar.fValues(:,2));            
            s.fValues(:,1) = (b.fValues(:,1) - bBar.fValues(:,1)).^2/bNorm1;
            s.fValues(:,2) = (b.fValues(:,2) - bBar.fValues(:,2)).^2/bNorm2;
            s.mesh = obj.mesh;
            s.order = 'P1';
            lRes = LagrangianFunction(s);

            %[iX,iY,~,~] = obj.computeIndex();
            %resX = res(iX,1);
            %resY = res(iY,1);
            %s.fValues = abs([resX,resY]);
            %s.mesh    = obj.mesh;
            %lRes = P1Function(s);
        end

        function nB = L2norm(obj,b)
            nB  = b'*(obj.massMatrix*b);
        end

        function [hRes,lRes] = evaluateAllResiduals(obj,rho,bBar,b)
            res = obj.evaluateResidual(rho,bBar,b);
            lRes = obj.evaluateLossResidual(res,bBar,b);            
            hRes = obj.evaluateHarmonicResidual(res,b);
        end

        function res = evaluateResidual(obj,rho,bBar,b)

           lambda = obj.computeInitalLambda();
           eta    = obj.computeInitialEta();
            %x = [b.fValues(:,1);b.fValues(:,2);lambda;eta];
           x = [b.fValues(:,1);b.fValues(:,2);eta];
           rhs = obj.computeRHS(rho,bBar);
           lhs = obj.computeLHS(rho,b);  
           res  = lhs*x - rhs;            
        end


        function [b,lambda] = solveProblem(obj,rho,bBar,bInitial)   
            b = bInitial;
            i = 1;
            hasNotConverged = true;
            lambda    = obj.computeInitalLambda();
            eta       = obj.computeInitialEta();
            etaOld    = eta;
            lambdaOld = lambda;
            rhs      = obj.computeRHS(rho,bBar);
            while hasNotConverged
             

                lhs = obj.computeLHS(rho,b);


                % gradient                            
                [iX,iY,iL,iE] = obj.computeIndex();
                %x = [b.fValues(:,1);b.fValues(:,2);lambda;eta];
                x = [b.fValues(:,1);b.fValues(:,2);eta];
                res  = lhs*x - rhs;
                resP = res([iX,iY]);
               % resD = res(iL);
                resE = res(iE);
           

                % Picard
                sol    = obj.solver.solve(lhs,rhs); 
                [iX,iY,iL,iE] = obj.computeIndex();
                bNew(:,1) = sol(iX,1);
                bNew(:,2) = sol(iY,1);
                lambda = sol(iL,1);
                eta    = sol(iE,1);

                
                
                err(i)  = norm(res)/norm(rhs);
                errP(i) = norm(resP)/norm(rhs);                
             %   errD(i) = norm(resD)/norm(rhs);
                errE(i) = norm(resE)/norm(rhs);
                incX(i) = norm(bNew(:)-b.fValues(:))/norm(bNew(:));
                incL(i) = norm(lambda-lambdaOld)/norm(lambdaOld);
                incE(i) = norm(eta-etaOld)/norm(etaOld);

                lambdaOld = lambda;
                etaOld = eta;
                errT(i) = max([err(i), incX(i), incL(i),incE(i) ]);
                hasNotConverged = errT(i) > 1e-3;%1e-13;

            

                 if mod(i,5) == 0
%                     obj.plotOrientation(b,2);                    
                     figure(99)%
                     %semilogy(1:i,([err;errP;errD;errE;incX;incL;incE]))                     
                     %legend('res','resP','resD','resE','incX','incL','incE')

                     semilogy(1:i,([err;errP;errE;incX;incL;incE]))
                     legend('res','resP','resE','incX','incL','incE')
%                     errN = obj.computeErrorNorm(b);
%                     obj.plotResiudalUnitBall(resE,errN)
                    figure(100)
                    b.plotArrowVector()
                    figure(101)
                    obj.plotSingulairties(b)

                    lRes = obj.evaluateLossResidual(res,bBar,b);
                    hRes = obj.evaluateHarmonicResidual(res,b);
                    hRes.plot
                    lRes.plot

                    quad = Quadrature.set(obj.mesh.type);            
                    quad.computeQuadrature('LINEAR');
                    xV = quad.posgp;

                    gB = b.evaluateGradient(xV);
                    gBV = gB.fValues;

           
                    s.fValues = abs(gBV([1:4],:)');
                    s.mesh    = obj.mesh;
                    s.order   = 'P0';
                    gb1 = LagrangianFunction(s);
                    gb1.plot
                    drawnow

          %           b = b.project('H1P1');
%
           %         b = obj.projectUnitBall(b.fValues);


                 end                
                                
                
                theta = 0.1;%0.5;
                bNew = obj.projectUnitBall(bNew);
                bNew = bNew.fValues;

                b = theta*bNew + (1-theta)*b.fValues ;    
                s.fValues = b;
                s.mesh = obj.mesh;
                s.order = 'P1';
                b = LagrangianFunction(s);

                
                i = i + 1;
            end



        end

        function a1 = createHalfOrientationVectorP1(obj,b1)
            bX = b1.fValues(:,1);
            bY = b1.fValues(:,2);
            beta   = atan2(bY,bX);
            alpha  = beta/2;
            aV(:,1) = cos(alpha);
            aV(:,2) = sin(alpha);
            s.fValues = aV;
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            a1 = LagrangianFunction(s);
        end

        function plotSingulairties(obj,b)
            a1   = obj.createHalfOrientationVectorP1(b);            
            s.mesh        = obj.mesh;
            s.orientation = a1;
            sC = SingularitiesComputer(s);
            sC.compute();
            sC.plot();
        end


        function plotHarmonicity(obj,resD)
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            resDT = zeros(size(x));
            iDOFs = obj.internalDOFs;
            resDT(iDOFs) = abs(resD);
            figure(59)
            clf
            trisurf(obj.mesh.connec,x,y,resDT)
            view(0,90)
            colorbar            
        end

        function [iX,iY,iL,iE] = computeIndex(obj)           
            npnod = obj.mesh.nnodes;
            iDOFs = obj.internalDOFs;
            iX = 1:npnod;
            iY = (npnod) + (1:npnod);
            iL = (2*npnod + 1):(2*npnod +length(iDOFs));
            %iE = (2*npnod + length(iDOFs) + 1):(2*npnod + length(iDOFs) + npnod);
            iE = (2*npnod +  1):(2*npnod + npnod);
        end

        function norm = computeNorm(obj,v)
         vx = v(:,1);
         vy = v(:,2);
         norm = sqrt(vx.^2 + vy.^2);
        end

        function tp = projectUnitBall(obj,t)
            uP = UnitBallProjector([]);
            tp = uP.project(t);
            s.fValues = tp;
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            tp = LagrangianFunction(s);
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
           obj.l1 = 100;%(100*obj.mesh.computeMeanCellSize).^2;                        
           obj.l2 = 100;           
        end

        function createInternalDOFs(obj)
           b     = obj.boundaryMesh;
           iDOFs = setdiff(1:obj.mesh.nnodes,b); 
           obj.internalDOFs = iDOFs;
        end
        
        function computeMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = LagrangianFunction.create(obj.mesh, 1, 'P1');
            s.trial = LagrangianFunction.create(obj.mesh, 1, 'P1');
            s.quadratureOrder = 'QUADRATICMASS';
            lhs = LHSintegrator.create(s);
            M = lhs.compute();           
            obj.massMatrix = M;
        end

        function computeStiffnessMatrix(obj)        
            s.test  = LagrangianFunction.create(obj.mesh, 1, 'P1');
            s.trial = LagrangianFunction.create(obj.mesh, 1, 'P1');
            s.mesh         = obj.mesh;
            s.type         = 'StiffnessMatrix';
            lhs = LHSintegrator.create(s);
            K = lhs.compute();
            obj.stiffnessMatrix = K;
        end

        function [CX,CY,DX,DY,EX,EY,Kxx,Kxy,Kyy,Mxx,Mxy,Myy,Mrho] = computeAdvectionMatrix(obj,rho,vH)
            s.fValues = zeros(obj.mesh.nnodes,1);
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            test = LagrangianFunction(s);
            
            s.mesh            = obj.mesh;
            s.type            = 'AdvectionMatrix';
            s.test            = test;
            s.trial           = test;
            s.b               = vH;
            s.rho             = rho;
            s.quadratureOrder = 'CUBIC';
            
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


       function  LHS = computeLHS(obj,rho,vH)
           l1 = obj.l1;
           l2 = obj.l2;
          
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
           Ex2 = diag(vH.fValues(:,1));
           Ey2 = diag(vH.fValues(:,2));
   %        lhs  = [Mrho+l2*(Kxx+Mxx)+K,Zb-l2*(Kxy+Mxy),(-Dx+Cx),Ex;...
   %                Zb-l2*(Kxy+Mxy),Mrho+l2*(Kyy+Myy)+K,(Dy-Cy) ,Ey;...
   %               (-Dx+Cx)',(Dy-Cy)',Z,Zbred';...
   %                Ex2',Ey2',Zbred,Zb];
           %lhs  = [M,Zb,(-Dx+Cx),Ex;Zb,M,(Dy-Cy),Ey;(-Dx+Cx)',(Dy-Cy)',Z,Zbred';Ex',Ey',Zbred,Zb];
      %     LHS = lhs;

           LHS  = [Mrho+l1*K+l2*(Kxx+Mxx),Zb-l2*(Kxy+Mxy),Ex;...
                   Zb-l2*(Kxy+Mxy),Mrho+l1*K+l2*(Kyy+Myy),Ey;...                  
                   Ex2',Ey2',Zb];
       end

       function Z = computeZeroMatrix(obj)
           iDOFs = obj.internalDOFs;
           Z     = zeros(length(iDOFs),length(iDOFs));           
       end

        function rhsV = computeRHS(obj,rho,bBar)
%             bBarV = squeeze(bBar.fValues);
% 
%     %         RHS not working for P0
%              bB = bBar.project('P0');
% % 
%              s.fValues = squeeze(bB.fValues(1,:,:));
%              s.mesh     = obj.mesh;
%              b1 = P0Function(s);
% % 
%              s.mesh = obj.mesh;
%              s.fValues = zeros(obj.mesh.nnodes,1);
%              test = P1Function(s);
%              s.mesh  = obj.mesh;
%              s.type = 'ShapeFunction';
%              rhs    = RHSintegrator.create(s);
%              rhs1   = rhs.compute(b1,test);
%              rhs2   = rhs.compute(b2,test);
%            
            rhoV = rho.fValues;
            l1 = obj.l1;

            %bBar1 = bBar.project('P1');
            w = 4*(1-rhoV).*rhoV;
            M = obj.massMatrix;            
            rhs1 = M*(bBar.fValues(:,1).*w);
            rhs2 = M*(bBar.fValues(:,2).*w);


            iDOFs = obj.internalDOFs;
            Z     = zeros(length(iDOFs),1);
            I     = ones(size(M,1),1);

%            rhs = [rhs1;rhs2;Z;M*I];
            %M2 = diag(sum(M));
           % rhs = [rhs1;rhs2;Z;M*I];
          
           
           % rhsV = [rhs1;rhs2;Z;I];
            rhsV = [rhs1;rhs2;I];


        end

      

    end
    
end