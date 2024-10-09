classdef EIFEMtesting < handle

    properties (Access = public)

    end
    properties (Access = private)
        meshDomain
        boundaryConditions
        bcApplier
        LHS
        RHS

    end


    properties (Access = private)
        nSubdomains
    end    

    methods (Access = public)



        function obj = EIFEMtesting()
            obj.init()

            close all            
            mR = obj.createReferenceMesh(); 
            bS  = mR.createBoundaryMesh();            
            [mD,iC,lG] = obj.createMeshDomain(mR);
            obj.meshDomain = mD;            
            [bC,dir] = obj.createBoundaryConditions(obj.meshDomain);
            obj.boundaryConditions = bC;
            obj.createBCapplier()
     
%             obj.createModelPreconditioning();
%             u = obj.solver2(LHS,RHS,refLHS);
            [LHS,RHS] = obj.createElasticProblem();
            s.LHS = LHS;
            s.RHS = RHS;
            s.meshDomain = obj.meshDomain;
            s.nSubdomains = obj.nSubdomains;
            s.interfaceConnec = iC;
            s.locGlobConnec   = lG;     
            s.nBoundaryNodes = bS{1}.mesh.nnodes;
            s.nReferenceNodes = mR.nnodes;
            s.coarseMesh      = obj.createCoarseMesh(mR);
            s.dir = dir;
            s.bcApplier = obj.bcApplier;
            gP = GeneralPreconditioner(s);


            LHSf = @(x) LHS*x;
            RHSf = RHS;      

            Usol = LHS\RHS;

            Mid = @(r) r;
            MidOrth = @(r,A,z) z+0.3*(r-A(z));

            Meifem = @(r) gP.solveEIFEM(r);
            s.LHS = LHS;
            P = PreconditionerILU(s);            
            Milu = @(r) P.apply(r);



            s.LHS = LHS;
            P = PreconditionerGaussSeidel(s);              
            MgaussSeidel = @(r) P.apply(r);

            s.LHS = LHS;
            P = PreconditionerJacobi(s);              
            MJacobi = @(r) P.apply(r);            

            MiluCG = @(r) gP.InexactCG(r,LHSf,MJacobi);            

         %  LHSf = @(x) P*LHS*x;            
            %  RHSf = P*RHS;
            tol = 1e-8;
            P = @(r) Mid(r); %obj.multiplePrec(r,Mid,Mid,LHSf);
            tic
            x0 = zeros(size(RHSf));
            [uCG,residualCG,errCG,errAnormCG] = PCG.solve(LHSf,RHSf,x0,P,tol,Usol);
            toc
            %[uCG,residualCG,errCG,errAnormCG] = RichardsonSolver.solve(LHSf,RHSf,x0,P,tol,0.1,Usol);
            
            
            M = MiluCG;%Milu_m;%Meifem; %Milu %Pm
            M2 = Meifem;
            M3 = MiluCG;
%             [uPCG,residualPCG,errPCG,errAnormPCG] = obj.solverTestEifem(LHSf,RHSf,Usol,M);
            tol = 1e-8;
            P = @(r) gP.multiplePrec(r,M,M2,M3,LHSf);  
%            P = Milu;
%              P = @(r) obj.additivePrec(r,Mid,Mmodal,LHSf);  
            tic
            x0 = zeros(size(RHSf));            
            [uPCG,residualPCG,errPCG,errAnormPCG] = PCG.solve(LHSf,RHSf,x0,P,tol,Usol);
            toc
            %[uCG,residualCG,errCG,errAnormCG] = RichardsonSolver.solve(LHSf,RHSf,x0,P,tol,0.1,Usol);

            figure
            plot(residualPCG,'linewidth',2)
            hold on
            plot(residualCG,'linewidth',2)
            set(gca, 'YScale', 'log')
            legend({'CG + ILU-EIFEM-ILU','CG'},'FontSize',12)
            xlabel('Iteration')
            ylabel('Residual')

             figure
            plot(errPCG,'linewidth',2)
            hold on
            plot(errCG,'linewidth',2)
            set(gca, 'YScale', 'log')
            legend('CG + EIFEM+ ILU(CG-90%-L2)','CG')
            xlabel('Iteration')
            ylabel('||error||_{L2}')

            figure
            plot(errAnormPCG,'linewidth',2)
            hold on
            plot(errAnormCG,'linewidth',2)
            set(gca, 'YScale', 'log')
            legend('CG + EIFEM+ ILU(CG-90%-L2)','CG')
            xlabel('Iteration')
            ylabel('Energy norm')

        end
 
    end

    methods (Access = private)

        function init(obj)
            obj.nSubdomains  = [5 1]; %nx ny
        end        

        function [mD,iC,lG] = createMeshDomain(obj,mR)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = mR;        
            m = MeshCreatorFromRVE(s);
            [mD,~,iC,~,lG] = m.create();
       end


       function mS = createReferenceMesh(obj)
        %    mS = obj.createStructuredMesh();
        %   mS = obj.createMeshFromGid();
           mS = obj.createEIFEMreferenceMesh();
       end       
       

       function mS = createMeshFromGid(obj)
           filename   = 'lattice_ex1';
           a.fileName = filename;
           femD       = FemDataContainer(a);
           mS         = femD.mesh;
       end

       function mS = createStructuredMesh(obj)

           % Generate coordinates
           x1 = linspace(0,1,5);
           x2 = linspace(0,1,5);
           % Create the grid
           [xv,yv] = meshgrid(x1,x2);
           % Triangulate the mesh to obtain coordinates and connectivities
           [F,coord] = mesh2tri(xv,yv,zeros(size(xv)),'x');

           s.coord    = coord(:,1:2);
           s.connec   = F;
           mS         = Mesh.create(s);
       end

       function mS = createEIFEMreferenceMesh(obj)
            filename = 'DEF_Q4porL_1.mat';
            load(filename);
            s.coord    = EIFEoper.MESH.COOR;
            s.coord(s.coord==min(s.coord)) = round(s.coord(s.coord==min(s.coord)));
            s.coord(s.coord==max(s.coord)) = round(s.coord(s.coord==max(s.coord)));
            s.connec   = EIFEoper.MESH.CN;
            mS         = Mesh.create(s);
       end

       function mCoarse = createCoarseMesh(obj,mR)
           s.nsubdomains   = obj.nSubdomains; %nx ny
           s.meshReference = obj.createReferenceCoarseMesh(mR);
           mRVECoarse      = MeshCreatorFromRVE(s);
           [mCoarse,~,~] = mRVECoarse.create();
       end



       function cMesh = createReferenceCoarseMesh(obj,mR)
           xmax = max(mR.coord(:,1));
           xmin = min(mR.coord(:,1));
           ymax = max(mR.coord(:,2));
           ymin = min(mR.coord(:,2));
           coord(1,1) = xmin;
           coord(1,2) = ymin;
           coord(2,1) = xmax;
           coord(2,2) = ymin;
           coord(3,1) = xmax;
           coord(3,2) = ymax;
           coord(4,1) = xmin;
           coord(4,2) = ymax;
           %             coord(1,1) = xmax;
           %             coord(1,2) = ymin;
           %             coord(2,1) = xmax;
           %             coord(2,2) = ymax;
           %             coord(3,1) = xmin;
           %             coord(3,2) = ymax;
           %             coord(4,1) = xmin;
           %             coord(4,2) = ymin;
           connec = [1 2 3 4];
           connec = [2 3 4 1];
           s.coord = coord;
           s.connec = connec;
           cMesh = Mesh.create(s);
       end

        function createBCapplier(obj)
            s.mesh                  = obj.meshDomain;
            s.boundaryConditions    = obj.boundaryConditions;
            obj.bcApplier           = BCApplier(s);
        end
       

        function material = createMaterial(obj,mesh)
            [young,poisson] = obj.computeElasticProperties(mesh);
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = mesh.ndim;
            s.young   = young;
            s.poisson = poisson;
            tensor    = Material.create(s);
            material  = tensor;
        end

        function [young,poisson] = computeElasticProperties(obj,mesh)
            E1  = 1;
            nu1 = 1/3;
            E   = AnalyticalFunction.create(@(x) E1*ones(size(squeeze(x(1,:,:)))),1,mesh);
            nu  = AnalyticalFunction.create(@(x) nu1*ones(size(squeeze(x(1,:,:)))),1,mesh);
            young   = E;
            poisson = nu;
        end        

        function [Dir,PL] = createRawBoundaryConditions(obj)
            isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
            isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
            isBottom  = @(coor) (abs(coor(:,2) - min(coor(:,2)))   < 1e-12);
            isTop  = @(coor) (abs(coor(:,2) - max(coor(:,2)))   < 1e-12);
            %             isMiddle = @(coor) (abs(coor(:,2) - max(coor(:,2)/2)) == 0);
            Dir{1}.domain    = @(coor) isLeft(coor)| isRight(coor) ;
            Dir{1}.direction = [1,2];
            Dir{1}.value     = 0;

            %             Dir{2}.domain    = @(coor) isRight(coor) ;
            %             Dir{2}.direction = [2];
            %             Dir{2}.value     = 0;

            PL.domain    = @(coor) isTop(coor);
            PL.direction = 2;
            PL.value     = -0.1;
        end

        function [bc,Dir,PL] = createBoundaryConditions(obj,mesh)
            [Dir,PL]  = obj.createRawBoundaryConditions();
            dirichletFun = [];
            for i = 1:numel(Dir)
                dir = DirichletCondition(obj.meshDomain, Dir{i});
                dirichletFun = [dirichletFun, dir];
            end

            pointload = PointLoad(mesh,PL);
            % need this because force applied in the face not in a point
            pointload.values        = pointload.values/size(pointload.dofs,1);
            fvalues                 = zeros(mesh.nnodes*mesh.ndim,1);
            fvalues(pointload.dofs) = pointload.values;
            fvalues                 = reshape(fvalues,mesh.ndim,[])';
            pointload.fun.fValues   = fvalues;

            s.pointloadFun = pointload;
            s.dirichletFun = dirichletFun;
            s.periodicFun  =[];
            s.mesh         = mesh;
            bc             = BoundaryConditions(s);
        end


        function [LHSr,RHSr] = createElasticProblem(obj)
            u = LagrangianFunction.create(obj.meshDomain,obj.meshDomain.ndim,'P1');
            material = obj.createMaterial(obj.meshDomain);
            [lhs,LHSr] = obj.computeStiffnessMatrix(obj.meshDomain,u,material);
            RHSr       = obj.computeForces(lhs,u);
        end


        function [LHS,LHSr] = computeStiffnessMatrix(obj,mesh,dispFun,mat)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = mesh;
            s.test     = dispFun;
            s.trial    = dispFun;
            s.material = mat;
            s.quadratureOrder = 2;
            lhs = LHSintegrator.create(s);
            LHS = lhs.compute();
            LHSr = obj.bcApplier.fullToReducedMatrixDirichlet(LHS);            
        end

        function RHS = computeForces(obj,stiffness,u)
            s.type      = 'Elastic';
            s.scale     = 'MACRO';
            s.dim.ndofs = u.nDofs;
            s.BC        = obj.boundaryConditions;
            s.mesh      = obj.meshDomain;
            RHSint      = RHSintegrator.create(s);
            rhs         = RHSint.compute();
            % Perhaps move it inside RHSint?
            R           = RHSint.computeReactions(stiffness);
            RHS = rhs+R;
            RHS = obj.bcApplier.fullToReducedVectorDirichlet(RHS);            
        end        


    end

end
