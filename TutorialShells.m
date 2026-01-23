classdef TutorialShells < handle


    properties (Access = private)
        mesh
        young
        area
        shear
        inertia
        uFun
        thetaFun
        wFun
        bcU,bcT,bcW
    end

    methods (Access = public)

        function obj = TutorialShells()
            obj.createMesh()
            obj.createMaterialProperties()
            obj.createSolutionField()
            
            obj.createBoundaryConditions()
            LHS = obj.createLHS();
            RHS = obj.createRHS();

            x = LHS\RHS;

            nTheta = length(obj.computeFreeDofs(obj.bcT));
            nU     = length(obj.computeFreeDofs(obj.bcU));
            nW     = length(obj.computeFreeDofs(obj.bcW));            
            uF = x(1:nU,1);
            tF = x((nU+1):(nU+nTheta),1);
            wF = x((nU+nTheta+1):(nU+nTheta+nW),1);

            dofFT = obj.computeFreeDofs(obj.bcT);
            dofFU = obj.computeFreeDofs(obj.bcU);
            dofFW = obj.computeFreeDofs(obj.bcW);

            uT = zeros(obj.uFun.nDofs,1);
            uT(dofFT,1) = uF; 
            uT = reshape(uT,[], obj.uFun.ndimf);
            obj.uFun.setFValues(uT);
            plot(obj.uFun)
        end

    end

    methods (Access = private)

        function createMesh(obj)
          obj.mesh = UnitTriangleMesh(50,50);
        end

        function createSolutionField(obj)
           obj.uFun     = LagrangianFunction.create(obj.mesh,2,'P1');
           obj.thetaFun = LagrangianFunction.create(obj.mesh,2,'P1');
           obj.wFun     = LagrangianFunction.create(obj.mesh,1,'P1');
        end

        function createMaterialProperties(obj)
          E = 3;
          obj.young = ConstantFunction.create(E,obj.mesh);
          obj.area = ConstantFunction.create(1,obj.mesh);
          obj.shear = ConstantFunction.create(1,obj.mesh);
          obj.inertia = ConstantFunction.create(1,obj.mesh);
        end

        function LHS = createLHS(obj)
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
            LHS = [Ku Zut Zuw; Zut' (Ktheta+Mtheta) Nthetaw; Zuw' Nthetaw' Kw];
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
            isRight = @(coor)  abs(coor(:,1)-xMax)< TOL;
            isTop   = @(coor)  abs(coor(:,2)-yMax)< TOL;
            isBotom = @(coor)  abs(coor(:,2)-yMin)< TOL;

            sDir{1}.domain    = @(coor) isLeft(coor) | isRight(coor) | isTop(coor) | isBotom(coor);
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





    end

end