classdef TutorialShellsDinamica < handle


    properties (Access = private)
        mesh
        young
        poisson
        density
        area
        shear
        inertia
        uFun
        thetaFun
        wFun
        bcU,bcT,bcW
    end

    methods (Access = public)

        function obj = TutorialShellsDinamica()
            obj.createMesh()
            obj.createMaterialProperties()
            obj.createSolutionField()
            
            obj.createBoundaryConditions()
            [LHS, mLHS] = obj.createLHS();
            RHS = obj.createRHS();

            [PHI,w2]=eigs(LHS,mLHS,10,'sm');
            sqrt(diag(w2))

            nTheta = length(obj.computeFreeDofs(obj.bcT));
            nU     = length(obj.computeFreeDofs(obj.bcU));
            nW     = length(obj.computeFreeDofs(obj.bcW)); 

            dofFT = obj.computeFreeDofs(obj.bcT);
            dofFU = obj.computeFreeDofs(obj.bcU);
            dofFW = obj.computeFreeDofs(obj.bcW);

            %for nf = 1:10
            %    wF = PHI((nU+nTheta+1):(nU+nTheta+nW),nf);
            %    wT = zeros(obj.wFun.nDofs,1);
            %    wT(dofFW,1) = wF; 
            %    wT = reshape(wT,[], obj.wFun.ndimf);
            %    obj.wFun.setFValues(wT);
            %    figure(nf)
            %    plot(obj.wFun);
            %end
            
            t  = linspace(0,500,500);
            %timeFunct = @(tau) heaviside(tau);
            
            timeFunct = @(tau) sin(0.5*tau);
            
            RHSnorm = PHI'*RHS;

            RHSt = RHSnorm*timeFunct(t);
            DM = PHI'*mLHS*PHI;
            DC = 0.01*eye(size(DM));
            DK = PHI'*LHS*PHI;

            q0    = zeros(size(DM,1),1);
            qdot0 = zeros(size(DM,1),1);


            [q, qdot, q2dot] = obj.NewmarkIntegration(DM,DC,DK,RHSt,t,q0,qdot0);

            sol = PHI*q;

            figure(1)
            xlim([0,t(end)]);
            hold on
            for tau=1:length(t)
                wF = sol((nU+nTheta+1):(nU+nTheta+nW),tau);
                wT = zeros(obj.wFun.nDofs,1);
                wT(dofFW,1) = wF; 
                wT = reshape(wT,[], obj.wFun.ndimf);
                obj.wFun.setFValues(wT);

                plot(t(tau),wT(3383),'ok');
                drawnow
                pause(0.05);

            end
            

        end

    end

    methods (Access = private)

        function createMesh(obj)
          %fullmesh = UnitTriangleMesh(50,50);
          fullmesh = TriangleMesh(18,10,60,60);
          ls = obj.computeWingLevelSet(fullmesh);
          sUm.backgroundMesh = fullmesh;
          sUm.boundaryMesh   = fullmesh.createBoundaryMesh;
          uMesh              = UnfittedMesh(sUm);
          uMesh.compute(ls);
          wingMesh = uMesh.createInnerMesh();
          obj.mesh = wingMesh;

          %obj.mesh = UnitTriangleMesh(50,50);

        end

        function ls = computeWingLevelSet(obj, mesh)
            gPar.type          = 'WingShape';
            gPar.xCoorCenter   = -0.01;
            gPar.yCoorCenter   = -0.01;
            gPar.chordRoot     = 7.3;
            gPar.chordTip      = 1.25;
            gPar.semiSpan      = 18.0;
            gPar.sweepDeg      = 25.0;
            g                  = GeometricalFunction(gPar);
            phiFun             = g.computeLevelSetFunction(mesh);
            lsWing             = phiFun.fValues;
            ls = lsWing;
        end

        function createSolutionField(obj)
           obj.uFun     = LagrangianFunction.create(obj.mesh,2,'P1');
           obj.thetaFun = LagrangianFunction.create(obj.mesh,2,'P1');
           obj.wFun     = LagrangianFunction.create(obj.mesh,1,'P1');
        end

        function createMaterialProperties(obj)
            E = 3;
            nu = 0.3;
            Density = 0.1;
            obj.young   = ConstantFunction.create(E, obj.mesh);
            obj.poisson = ConstantFunction.create(nu, obj.mesh);
            obj.density = ConstantFunction.create(Density, obj.mesh);
            obj.area    = ConstantFunction.create(1,obj.mesh);
            obj.shear   = ConstantFunction.create(1,obj.mesh);
            obj.inertia = ConstantFunction.create(1/12,obj.mesh);

        end

        function [LHS, mLHS] = createLHS(obj)
            E = obj.young;
            A = obj.area;
            f = @(u,v) E.*A.*DDP(SymGrad(v),SymGrad(u));
            Ku = IntegrateLHS(f,obj.uFun,obj.uFun,obj.mesh,'Domain',2);
            Ku = obj.reduceMatrix(Ku,obj.bcU,obj.bcU);

            E = obj.young;
            I = obj.inertia;
            f = @(u,v) E.*I.*DDP(SymGrad(v),SymGrad(u));
            Ktheta = IntegrateLHS(f,obj.thetaFun,obj.thetaFun,obj.mesh,'Domain',2);
            Ktheta = obj.reduceMatrix(Ktheta,obj.bcT,obj.bcT);


            A = obj.area;
            G = obj.shear;
            f = @(u,v) G.*A.*DP(v,u);
            Mtheta = IntegrateLHS(f,obj.thetaFun,obj.thetaFun,obj.mesh,'Domain',2);
            Mtheta = obj.reduceMatrix(Mtheta,obj.bcT,obj.bcT);


            A = obj.area;
            G = obj.shear;
            f = @(u,v) G.*A.*DP(v,Grad(u));
            Nthetaw = IntegrateLHS(f,obj.thetaFun,obj.wFun,obj.mesh,'Domain',2);            
            Nthetaw = obj.reduceMatrix(Nthetaw,obj.bcT,obj.bcW);

        
            A = obj.area;
            G = obj.shear;
            f = @(u,v) A.*G.*DP(Grad(v),Grad(u));
            Kw = IntegrateLHS(f,obj.wFun,obj.wFun,obj.mesh,'Domain',2);  
            Kw = obj.reduceMatrix(Kw,obj.bcW,obj.bcW);

            nTheta = length(obj.computeFreeDofs(obj.bcT));
            nU     = length(obj.computeFreeDofs(obj.bcU));
            nW     = length(obj.computeFreeDofs(obj.bcW));

            Zut = zeros(nU,nTheta);
            Zuw = zeros(nU,nW);
            LHS = [Ku Zut Zuw; 
                   Zut' (Ktheta+Mtheta) Nthetaw;
                   Zuw' Nthetaw' Kw];


            rho = obj.density;
            f = @(u,v) rho.*DP(u,v);
            Mu = IntegrateLHS(f,obj.uFun,obj.uFun,obj.mesh,'Domain',2);
            Mu = obj.reduceMatrix(Mu,obj.bcU,obj.bcU);

            rho = obj.density;
            f = @(u,v) rho.*DP(u,v);
            Mw = IntegrateLHS(f,obj.wFun,obj.wFun,obj.mesh,'Domain',2);
            Mw = obj.reduceMatrix(Mw,obj.bcW,obj.bcW);

            rho = obj.density;
            f = @(u,v) rho.*DP(u,v);
            Mtheta = IntegrateLHS(f,obj.uFun,obj.uFun,obj.mesh,'Domain',2);
            Mtheta = obj.reduceMatrix(Mtheta,obj.bcT,obj.bcT);


            Ztw = zeros(nTheta,nW);
            mLHS = [Mu, Zut, Zuw;
                   Zut' Mtheta Ztw;
                   Zuw' Ztw' Mw];

        end

        function RHS = createRHS(obj)
            p = ConstantFunction.create([0 0],obj.mesh);
            m = ConstantFunction.create([0 0],obj.mesh);
            q = ConstantFunction.create(1,obj.mesh);

            fu = @(v) DP(p,v);
            RHSu = IntegrateRHS(fu,obj.uFun,obj.mesh,'Domain',2);
            RHSu = obj.reduceVector(RHSu,obj.bcU);

            ftheta   = @(v) DP(m,v);
            RHStheta = IntegrateRHS(ftheta,obj.thetaFun,obj.mesh,'Domain',2);
            RHStheta = obj.reduceVector(RHStheta,obj.bcT);

            fw = @(v) q.*v;
            RHSw = IntegrateRHS(fw,obj.wFun,obj.mesh,'Domain',2);
            RHSw = obj.reduceVector(RHSw,obj.bcW);
            

            RHS = [RHSu;RHStheta;RHSw];
        end

        function createBoundaryConditions(obj)
            obj.bcU = obj.createGeneralBoundaryConditions([1 2]);
            obj.bcT = obj.createGeneralBoundaryConditions([1 2]);
            obj.bcW = obj.createGeneralBoundaryConditions([1]);
        end

        function bc = createGeneralBoundaryConditions(obj,direct)
            TOL = 1e-12;
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            xMin    = min(obj.mesh.coord(:,1));
            yMin    = min(obj.mesh.coord(:,2));
            isLeft  = @(coor)  abs(coor(:,1)-xMin)< TOL;
            %isRight = @(coor)  abs(coor(:,1)-xMax)< TOL;
            isTop   = @(coor)  abs(coor(:,2)-yMax)< TOL;
            isBotom = @(coor)  abs(coor(:,2)-yMin)< TOL;

            sDir{1}.domain    = @(coor) isLeft(coor);
            %sDir{1}.domain    = @(coor) isLeft(coor) | isRight(coor) | isTop(coor) | isBotom(coor);
            
            sDir{1}.direction = direct;
            sDir{1}.value     = 0;                    
            sDir{1}.ndim = length(direct);


            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;

            s.periodicFun  = [];
            s.pointloadFun    = [];
            s.mesh = obj.mesh;
            bc = BoundaryConditions(s);            
        end

        function fD = computeFreeDofs(obj,bC)
            dofs = 1:bC.dirichletFun.nDofs;
            fD   = setdiff(dofs,bC.dirichlet_dofs); 
        end

        function LHSred = reduceMatrix(obj,LHS,bcV,bcU)
            fdofV = obj.computeFreeDofs(bcV);
            fdofU = obj.computeFreeDofs(bcU);
            LHSred = LHS(fdofV,fdofU);
        end

        function RHSred = reduceVector(obj,RHS,bc)
            fdofV = obj.computeFreeDofs(bc);
            RHSred = RHS(fdofV,1);
        end

        function [q,qdot,q2dot] = NewmarkIntegration(obj,M,C,K,Ft,t,q0,qdot0)
            q       = zeros(size(K,1),length(t));
            qdot    = zeros(size(K,1),length(t));
            q2dot = zeros(size(K,1),length(t));

            alpha = 0.5;
            gamma = 0.5;
            dt = t(2)-t(1);

            a1 = (1-alpha)*dt;
            a2 = alpha*dt;
            a3 = 2/(gamma*dt^2);
            a4 = a3*dt;
            a5 = (1-gamma)/gamma;

            b0 = a3;
            b1 = a4;
            b2 = a5;
            b3 = a2*a3;
            b4 = a2*a4-1;
            b5 = a2*a5-a1;

            A1 = inv(b0*M+b3*C+K);
            A2 = A1*M;
            A3 = A1*C;

            q(:,1)=q0;
            qdot(:,1)=qdot0;
            q2dot(:,1)=inv(M)*Ft(:,1);

            for n = 2 : length(t)-1
                q(:,n) = A1*Ft(:,n)+...
                        A2*(b0*q(:,n-1)+b1*qdot(:,n-1)+b2*q2dot(:,n-1))+...
                        A3*(b3*q(:,n-1)+b4*qdot(:,n-1)+b5*q2dot(:,n-1));

                q2dot(:,n) = a3*(q(:,n)-q(:,n-1))-a4*qdot(:,n-1)-a5*q2dot(:,n-1);
                qdot(:,n) = qdot(:, n-1) +a1*q2dot(:,n-1)+a2*q2dot(:,n);

            end





        end





    end

end