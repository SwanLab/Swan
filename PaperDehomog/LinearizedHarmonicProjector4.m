classdef LinearizedHarmonicProjector4 < handle

    properties (Access = private)
        eta
        internalDOFs
        massMatrixSS
        massMatrixGG
        massMatrixBB
        stiffnessMatrixBB
        fB
        fBV
        fS
        perimeter
        filter
    end

    properties (Access = private)
        mesh
        boundaryNodes
        density
    end

    methods (Access = public)

        function obj = LinearizedHarmonicProjector4(cParams)
            obj.init(cParams);
            obj.fB = LagrangianFunction.create(obj.mesh, 1, 'P1');
            obj.fBV = LagrangianFunction.create(obj.mesh, 2, 'P1');
            obj.fS = LagrangianFunction.create(obj.mesh, 1, 'P1');
            obj.createInternalDOFs();                      
            obj.createFilter();
            obj.eta = (0*obj.mesh.computeMeanCellSize)^2;
            rhoEps = 1e-3;
            rhoMin = min(obj.density.fValues) - rhoEps; 
            rhoMax = max(obj.density.fValues) + rhoEps; 
            %obj.perimeter = 4*(obj.density-rhoMin).*(rhoMax-obj.density);
            obj.perimeter = ConstantFunction.create(0.5,obj.mesh);%
            obj.massMatrixBB      = IntegrateLHS(@(u,v) DP(v,obj.perimeter.*u),obj.fB,obj.fB,obj.mesh,'Domain',4);
            obj.massMatrixSS      = IntegrateLHS(@(u,v) DP(v,u),obj.fS,obj.fS,obj.mesh,'Domain',4);
            obj.stiffnessMatrixBB = IntegrateLHS(@(u,v) DP(Grad(v),(1-obj.perimeter).*Grad(u)),obj.fB,obj.fB,obj.mesh,'Domain');            
        end


        function b = solveProblem(obj,bBar,b)
            b    = project(b,obj.fB.order);
            %bBar = project(bBar,obj.fB.order);
            RHS = obj.computeRHS(bBar);
            LHS = obj.computeLHS(b);
            nInt = size(obj.internalDOFs,2);
            x = [b.fValues(:);zeros(nInt,1)];
            res = norm(LHS*x - RHS)/norm(x);
            [resL,resH,resB,resG] = obj.evaluateResidualNorms(bBar,b);
            i = 1;
            thetaH = 0.01;
            thetaR = 1;
            thetaB = 1;
            

            nSing = 100;

            while res(i) > 1e-12 && i<=2
                
                %iter Harmonic 
                x   = LHS\RHS;
                bNew = obj.createVectorFromSolution(full(x));
                b    = obj.relaxationInSphere(bNew,b,thetaH);
               
              
                %iter Projection UnitBall
                bNew   = obj.projectInUnitBall(b);
                b = obj.relaxationInSphere(bNew,b,thetaB);                


               % sCf   = obj.createSingularities(b);
               % nSing = sum(sCf.fValues);
                %iter Filter
                if nSing >= 1
                bNew = obj.filterVector(b);
                b    = obj.relaxationInSphere(bNew,b,thetaR);
                end                
                
                %iter Projection UnitBall
                bNew   = obj.projectInUnitBall(b);
                b = obj.relaxationInSphere(bNew,b,thetaB); 

                sCf   = obj.createSingularities(b);
                nSing = sum(sCf.fValues);
               
                LHS = obj.computeLHS(b);
                RHS = obj.computeRHS(bBar);

                i   = i+1;
                res(i) = norm(LHS*x - RHS)/norm(x);
                [resL(i),resH(i),resB(i),resG(i)] = obj.evaluateResidualNorms(bBar,b);

              
              %  disp(['iter ', num2str(i),'  resL = ', num2str(resL(i)),'  resH = ', num2str(resH(i)),'  resB = ', num2str(resB(i)), '  resG = ', num2str(resG(i))])
                fprintf('iter %3d | nSing = %3d |  res = %10.4e | resL = %10.4e | resH = %10.4e | resB = %10.4e | resG = %10.4e\n', ...
                i, nSing, res(i), resL(i), resH(i), resB(i), resG(i));
                close all             


                 plotVector(b)
                 %obj.plotSingularities(b)
                 fig = figure(1);  set(fig, 'Units', 'normalized', 'OuterPosition', [0.5 0 0.5 1]);
                 %drawnow  
                  
                 %obj.plotgif(i-1,fig);

                
                 

            end
            plotVector(obj.projectInUnitBall(b));
            figure()
            plot(1:i,log10([res; resL; resH; resB; resG]))
            legend('LHS*x-RHS','resDistance','resHarmonic','resUnitBall','resGradient')

            % Mbb  = obj.massMatrixBB;
            % Kbb  = obj.stiffnessMatrixBB;            
            % A    = Mbb + obj.eta*Kbb;
            % rhsB = IntegrateRHS(@(v) DP(v,obj.perimeter.*bBar),bBar,obj.mesh,'Domain',4);   
            % Z  = sparse(obj.fB.nDofs,obj.fB.nDofs);            
            % xV    = [A Z; Z' A]\rhsB;
            % s.fValues = reshape(full(xV),[],2);
            % s.mesh    = obj.mesh;
            % s.order   = obj.fB.order;
            % b = LagrangianFunction(s);
            b = project(b,'P1');
        end

        function [resLnorm,resHnorm,resBnorm,resGnorm] = evaluateResidualNorms(obj,bBar,b)
            [resL,resH,resB,resG] = obj.evaluateAllResiduals(bBar,b);
            resLnorm = Norm(resL,'L2');
            resHnorm = Norm(resH,'L2');
            resBnorm = Norm(resB,'L2');
            resGnorm = Norm(resG,'L2');
        end

        function [resL,resH,resB,resG] = evaluateAllResiduals(obj,bBar,b)
            resL = project(DP(b-bBar,obj.perimeter.*(b-bBar)),'P1');
            resH = obj.evaluateHarmonicResidual(b);                
            resB = project(norm(b)- 1,'P1');
            resG = project(DDP(Grad(b),(1-obj.perimeter).*Grad(b)),'P1');
        end

        function resH = evaluateHarmonicResidual(obj,b)
            bs = b.getVectorFields();            
            f = (-Grad(bs{1}).*bs{2}+Grad(bs{2}).*bs{1});
            rhsV = IntegrateRHS(@(v) DP(Grad(v),f),obj.fS,obj.mesh,'Domain',4);
            rhsV(obj.boundaryNodes) = 0;
            Mss = obj.massMatrixSS;
            Kss = IntegrateLHS(@(u,v) DP(Grad(v),Grad(u)),obj.fS,obj.fS,obj.mesh,'Domain');
            eps = (1*obj.mesh.computeMeanCellSize)^2;
            LHS = Mss+eps*Kss;
            hf = LHS\rhsV;
            resH = obj.createFunction(full(hf),obj.fS.order);
        end

               


    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh             = cParams.mesh;
            obj.boundaryNodes    = cParams.boundaryMesh;
            obj.density          = cParams.density;            
        end
  
        function createInternalDOFs(obj)
            bNodes = obj.boundaryNodes;
            iDOFs  = setdiff(1:obj.fS.nDofs,bNodes);
            obj.internalDOFs = iDOFs;
        end

        function f = createFunction(obj,fV,order)
            s.fValues = fV;
            s.mesh    = obj.mesh;
            s.order   = order;
            f = LagrangianFunction(s);
        end

        function createFilter(obj)
            s.filterType   = 'PDE';
            s.mesh         = obj.mesh;
            s.trial        = obj.fB;
            f              = Filter.create(s);
            obj.filter     = f;            
        end

        function LHS = computeLHS(obj,b)
            Mbb  = obj.massMatrixBB;
            Kbb  = obj.stiffnessMatrixBB;
            nInt = size(obj.internalDOFs,2);
            bs = b.getVectorFields();
            normB  = (norm(b) - 1);
            MnormB = IntegrateLHS(@(u,v) DP(v,2*normB.*u),obj.fB,obj.fB,obj.mesh,'Domain',4);
            Kb1 = IntegrateLHS(@(u,v) DP(Grad(v),bs{1}.*Grad(u)),obj.fB,obj.fS,obj.mesh,'Domain',4);
            Kb2 = IntegrateLHS(@(u,v) DP(Grad(v),bs{2}.*Grad(u)),obj.fB,obj.fS,obj.mesh,'Domain',4);
            
            %Nb1 = computeAdvection(bs{1},obj.fB,obj.fB,obj.mesh,4);
            %Nb2 = computeAdvection(bs{2},obj.fB,obj.fB,obj.mesh,4);

            %v = LagrangianFunction.create(obj.mesh, 2, 'P1');

            %Nb1t = IntegrateLHS(@(u,v) u.*DP(Grad(v),Grad(bs{1})),obj.fB,obj.fS,obj.mesh,'Domain',4)';            
            %Nb2t = IntegrateLHS(@(u,v) u.*DP(Grad(v),Grad(bs{2})),obj.fB,obj.fS,obj.mesh,'Domain',4)'; 

            Nb1 = IntegrateLHS(@(u,v) DP(Grad(u),Grad(bs{1})).*v,obj.fB,obj.fS,obj.mesh,'Domain',4);
            Nb2 = IntegrateLHS(@(u,v) DP(Grad(u),Grad(bs{2})).*v,obj.fB,obj.fS,obj.mesh,'Domain',4);

            iDOFs = obj.internalDOFs;
            Kb1 = Kb1(:,iDOFs);
            Kb2 = Kb2(:,iDOFs);
            Nb1 = Nb1(:,iDOFs);
            Nb2 = Nb2(:,iDOFs);
            Z  = sparse(obj.fB.nDofs,obj.fB.nDofs);
            Zh = sparse(nInt,nInt);
            eps = (0*obj.mesh.computeMeanCellSize)^2;  
            A  = Mbb + obj.eta*Kbb + eps*MnormB;
            LHS = [A          ,          Z, (-Kb2+Nb2);...
                   Z'         ,          A,  (Kb1-Nb1);...
                   (-Kb2+Nb2)', (Kb1-Nb1)',         Zh];
        end

        function RHS = computeRHS(obj,bBar)
            rhsB = IntegrateRHS(@(v) DP(v,obj.perimeter.*bBar),obj.fBV,obj.mesh,'Domain',3);
            rhsB = reshape(rhsB,2,[])';
            rhsB = rhsB(:);
            rhsH = zeros(size(obj.internalDOFs,2),1);
            RHS  = [rhsB;rhsH];
        end

        function b = createVectorFromSolution(obj,x)
            nB = obj.fB.nDofs;
            bV = x(1:2*nB);
            s.fValues = reshape(bV,[],2);
            s.mesh    = obj.mesh;
            s.order   = obj.fB.order;
            b = LagrangianFunction(s);
        end        

        function bNew = filterVector(obj,b)
            %eps = obj.eta;
            eps = (20*obj.mesh.computeMeanCellSize)^2;
            obj.filter.updateEpsilon(eps);
            bVector = b.getVectorFields;
            bNew1 = obj.filter.compute(bVector{1},2);
            bNew2 = obj.filter.compute(bVector{2},2);
            bNew(:,1) = bNew1.fValues;
            bNew(:,2) = bNew2.fValues;
            s.fValues = reshape(bNew,[],2);
            s.mesh    = obj.mesh;
            s.order   = obj.fB.order;
            bNew = LagrangianFunction(s);
        end

        function sCf = createSingularities(obj,b)
            a1   = obj.createHalfOrientationVectorP1(b);
            s.mesh        = obj.mesh;
            s.orientation = a1;
            sC = SingularitiesComputer(s);
            sCf = sC.compute();
        end

        function plotSingularities(obj,b)
            sCf = obj.createSingularities(b);
            plot(sCf);            
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

    end

    methods (Access = private, Static)

        function x = relaxationInSphere(xNew,x,theta)
            phiG = ScalarProduct(xNew,x,'L2');
            w    = max(acos(phiG),1e-14);
            x = (sin((1-theta)*w)/sin(w)).*x + (sin(theta*w)/sin(w)).*xNew;
            %x = (1-theta).*x + theta.*xNew;
        end                 

        function vNF = projectInUnitBall(vF)
            v    = vF.fValues;
            norm = vecnorm(v,2,2);
            v    = v./norm;
            vNF = LagrangianFunction.create(vF.mesh,vF.ndimf,vF.order);
            vNF.setFValues(v);
        end


        function plotgif(i,fig)
            filename = 'anim.gif';
            frame = getframe(fig);
            [A,map] = rgb2ind(frame.cdata,256);
            if i==1
                imwrite(A,map,filename,"gif","LoopCount",0,"DelayTime",0.1);
            else
                imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",0.1);
            end
        end        

    end    

end