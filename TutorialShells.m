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
        bcU,bcT,bcW,bcq
        lhs,RHSS
        solverType
        type, values, fun
    end

    methods (Access = public)

        function obj = TutorialShells()
            clc; close all; 
            obj.createMesh()
            obj.createMaterialProperties()
            obj.createSolutionField()

            obj.solverType = 'REDUCED';
            
            obj.createBoundaryConditions()
            LHS = obj.createLHS();
            obj.lhs = LHS;
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
            

            wT = zeros(obj.wFun.nDofs,1);
            wT(dofFW,1) = wF; 
            wT = reshape(wT,[], obj.wFun.ndimf);
            obj.wFun.setFValues(wT);
            
            % wT = zeros(obj.wFun.nDofs,1);    CAMBIAR PARA THETA
            % wT(dofFW,1) = wF; 
            % wT = reshape(wT,[], obj.wFun.ndimf);
            % obj.wFun.setFValues(wT);
            
            plot(obj.wFun)
            plot(obj.uFun)
            plot(obj.thetaFun)
            
            obj.wFun.print('wfun print','Paraview')
            obj.uFun.print('ufun print','Paraview')
            obj.thetaFun.print('thetafun print','Paraview')
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
            q = ConstantFunction.create([0],obj.mesh);   

            fu = @(v) DP(p,v);
            RHSu = IntegrateRHS(fu,obj.uFun,obj.mesh,'Domain',2);
            RHSu = obj.reduceVector(RHSu,obj.bcU);

            ftheta   = @(v) DP(m,v);
            RHStheta = IntegrateRHS(ftheta,obj.thetaFun,obj.mesh,'Domain',2);
            RHStheta = obj.reduceVector(RHStheta,obj.bcT);

            fw = @(v) q.*v;
            RHSw = IntegrateRHS(fw,obj.wFun,obj.mesh,'Domain',2);
            RHSw = obj.reduceVector(RHSw,obj.bcW);

            obj.computeForces();
            RHSq = obj.RHSS;
            RHSq = obj.reduceVector(RHSq,obj.bcW);
            

            RHS = [RHSu;RHStheta;RHSq];
        end


        function computeForces(obj)
            bc  = obj.bcq;
            t   = bc.tractionFun;
            rhs = zeros(obj.wFun.nDofs,1);
            if ~isempty(t)
                for i = 1:numel(t)
                    rhsi = t(i).computeRHS(obj.wFun);
                    rhs  = rhs + rhsi;
                end
            end
            if strcmp(obj.solverType,'REDUCED')
                bc      = obj.bcq;
                dirich  = bc.dirichlet_dofs;
                dirichV = bc.dirichlet_vals;
                if ~isempty(dirich)
                    R = -obj.lhs(obj.wFun.nDofs(:),dirich)*dirichV; %??????
                else
                    R = zeros(sum(obj.wFun.nDofs(:)),1);
                end
                rhs = rhs+R;
            end
            obj.RHSS = rhs;
        end

        function createBoundaryConditions(obj)            
            obj.bcU = obj.createGeneralBoundaryConditions([1 2]);
            obj.bcT = obj.createGeneralBoundaryConditions([1 2]);
            obj.bcW = obj.createGeneralBoundaryConditions([1]);
            obj.bcq = obj.bcW;
            
        end


        % ============================================================================================================================================================
        % ============================================================================================================================================================

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


            % MODIFICATIONS: Embedded beam on one side ant with punctual load on the other  

            sDir{1}.domain    = @(coor) isLeft(coor);
            sDir{1}.direction = direct;
            sDir{1}.value     = 0; 
            sDir{1}.ndim      = length(direct);


            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;

            
            %% Apply point load to a single node on the right boundary
            
            % % Find right boundary nodes
            % rightNodes = find(isRight(obj.mesh.coord));
            % if isempty(rightNodes)
            %     error('No nodes found on the right boundary.');
            % end
            % % Choose one node: here pick the middle one (can change as needed)
            % idx = ceil(numel(rightNodes)/2);
            % singleNode = rightNodes(idx);
            % 
            % % Define point load structure for that single node
            % sPL{1}.domain    = @(coor) (1:size(coor,1))'==singleNode;
            % sPL{1}.direction = 2;
            % sPL{1}.value     = 1;

            %%

            sPL{1}.domain    = @(coor) isRight(coor);
            sPL{1}.direction = 2;
            sPL{1}.value     = 1;

            %%

            pointloadFun = [];
            for i = 1:numel(sPL)
                pl = TractionLoad(obj.mesh, sPL{i}, 'DIRAC');
                pointloadFun = [pointloadFun, pl];

            end
            s.pointloadFun = pointloadFun;


            s.periodicFun  = [];
            s.mesh = obj.mesh;
            bc = BoundaryConditions(s);            
        end


        % ============================================================================================================================================================
        % ============================================================================================================================================================

        function obj = TractionLoad2(mesh,s,type)
         
            obj.type = type;
            switch type
                case 'DIRAC'
                    pl_dofs = s.domain(mesh.coord);
                    vals = zeros(length(pl_dofs),1);
                    vals(pl_dofs,1) = s.value;
                    obj.values = reshape(vals',[],1);
                case 'FUNCTION'
                    dom     = s.domain;
                    f       = s.fun;
                    neuFun  = AnalyticalFunction.create(dom,f.mesh);
                    obj.fun = f.*neuFun;
                    obj.mesh = mesh;
            end
        end


        % ============================================================================================================================================================
        % ============================================================================================================================================================
        

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