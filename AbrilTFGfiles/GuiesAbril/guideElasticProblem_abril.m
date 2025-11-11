classdef guideElasticProblem_abril < handle

    % This code is used to learn how to compute the Elastic Problem for a
    % the coarse space. The mesh and geometry of the problem consists of a
    % square with a circular hole of radius r inputted to the program. 
    % There are slight variations from the original code, since it does not
    % take into account the scaling problem between the cases of different
    % r

    properties (Access = public)
        stiffness
        strain
        stress
        dLambda
        bcApplier
    end
    
    properties (Access = private)
        radius
        nodeDirection
        physicalProblem
    
        boundaryConditions
        problemSolver
        forces
        uFun
        strainFun
    
    end
    
    properties  (Access = protected)
        mesh
        material
        displacementFun
        boundaryMesh
        boundaryMeshJoined
        localGlobalConnecBd
    end

    methods (Access = public)

        function [obj,u,L]=guideElasticProblem_abril(r)
            obj.init(r, 1)
            obj.createMesh();
            [u, L] = obj.doElasticProblemHere();
            mesh = obj.mesh;

            z.mesh      = obj.mesh;
            z.order     = 'P1';

            for i=1:8
               z.fValues   = reshape(u(:,i),[obj.mesh.ndim,obj.mesh.nnodes])';
               uFeFun = LagrangianFunction(z);%
               fileName = ['trialGuide' num2str(i)];
               uFeFun.print(fileName,'Paraview');
            end
        end


    end

    methods (Access=private)
        function init (obj,r,i)
            obj.radius=r;
            obj.nodeDirection = i;
        end

        function createMesh(obj)
            bgMesh=obj.createReferenceMesh();  % Crea la background mesh
            lvSet=obj.createlvSetFunction(bgMesh);  % Crea la levelSet per fer el forat
            uMesh=obj.computeUnfittedMesh(bgMesh,lvSet);  %Fa el forat
            obj.mesh = uMesh.createInnerMesh();    % guarda com a mesh el conjunt 

            obj.boundaryMesh = obj.mesh.createBoundaryMesh();   %crea el boundary de la mesh
            [obj.boundaryMeshJoined, obj.localGlobalConnecBd] = obj.mesh.createSingleBoundaryMesh(); %conectivitats
        end


            function bgMesh=createReferenceMesh(~)
                % Generate coordinates and discretization
                x1 = linspace(-1,1,50);
                x2 = linspace(-1,1,50);
            
                % Create the grid
                [xv,yv] = meshgrid(x1,x2);
            
                % Triangulate the mesh to obtain coordinates and connectivities
                [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
                s.coord  = V(:,1:2);
                s.connec = F;
                bgMesh = Mesh.create(s);
            end
    
            function levelSet=createlvSetFunction(obj,bgMesh)
                sLS.type        = 'CircleInclusion';
                sLS.xCoorCenter = 0;    %coordenada y del centre
                sLS.yCoorCenter = 0;  %coordenada y del centre
                sLS.radius      = obj.radius;
                g               = GeometricalFunction(sLS);
                lsFun           = g.computeLevelSetFunction(bgMesh);
                levelSet        = lsFun.fValues;
            end
    
            function uMesh=computeUnfittedMesh(~,bgMesh,lvSet)
                sUm.backgroundMesh  = bgMesh;          % Defineix la background mesh
                sUm.boundaryMesh    = bgMesh.createBoundaryMesh();   %calcula el boundary mesh corresponent
                uMesh               = UnfittedMesh(sUm);   % background+boundary
                uMesh.compute(lvSet);     % inclusio radi
            end

        function [u,L]=doElasticProblemHere(obj)
            s=obj.createFEMContainer();        % Inicialitza les dades pel solver
            obj.createDisplacementFunHere();   % Crea la funcion de FE pels desplaçaments
           % obj.createBCApplyerHere(s);        % Crea les BC --> dubte de si es necessari
            obj.createSolverHere(s)            
            obj.computeStiffnessMatrixHere();  % crea la matriu K
           % obj.computeForcesHere(s);         % crea el vector F --> Aqui
                                               % no es necessari pq ell el calcula diferent

             c = obj.computeCmatP1();          % crea la matriu c
             rdir = obj.RHSdir();              % crea els vectors u_d
            [u, L]  = obj.computeDisplacementHere(c, rdir);
        end


            function s=createFEMContainer(obj)   %Aqui es defineixen les propietats necessaries pel FEM  sobre el solver
                s.mesh     = obj.mesh;  % malla a resoldre
                s.scale    = 'MACRO';   
                obj.createMaterial();   % Defineix les propietats del material, si es isotropic, elastic, etc
                s.material = obj.material;
                s.dim      = '2D';          % Dimensions del problema
               % s.boundaryConditions = obj.createBoundaryConditions(); % Revisar perque no crec que sigui necessari
                s.interpolationType  = 'LINEAR';
                s.solverType         = 'MONOLITHIC';
                s.solverMode         = 'DISP';
                s.solverCase         = 'DIRECT';
            end

                function createMaterial(obj,s)
                    [young,poisson] = obj.computeElasticProperties();
                    s.type       = 'ISOTROPIC';
                    s.ptype      = 'ELASTIC';
                    s.ndim       = obj.mesh.ndim;
                    s.young      = young;
                    s.poisson    = poisson;
                    obj.material = Material.create(s);            
                end

                %Aquesta funcio es per crear els 2 materials per tenir el mateix
                %scaling en tots els casos de radi, com en aquest exemple guia
                %nomes es fa per un radi no fa falta i es fara el exemple del
                %tutorial estandard
                function [young,poisson] = computeElasticProperties(obj)
                    E1  = 1;
%                   E2 = E1/1000;
                    nu = 1/3;
%                    x0=0;
%                    y0=0;
%                     young   = ConstantFunction.create(E,mesh);
%                     poisson = ConstantFunction.create(nu,mesh);
%                    f   = @(x) (sqrt((x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2)<obj.radius)*E2 + ...
%                                (sqrt((x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2)>=obj.radius)*E1 ; 
%                                              x(2,:,:).*0 ];
%                    young   = AnalyticalFunction.create(f,1,mesh);
                     young   = ConstantFunction.create(E1,obj.mesh);  
                     poisson = ConstantFunction.create(nu,obj.mesh);  
        
                end

                function bc=createBoundaryConditions(obj)
                    v                    = zeros(8,1);
                    v(obj.nodeDirection) = 1;
                    nRes                 = [1 1 2 2 3 3 4 4]*v;
                    assignMatrix         = [2 1 0 0 0 0 0 0
                                            0 0 2 1 0 0 0 0
                                            0 0 0 0 2 1 0 0
                                            0 0 0 0 0 0 2 1
                                            1 2 1 2 1 2 1 2];
        
                    vSimp       = assignMatrix*v;
                    dirs        = cell(5,1);
                    [dirs{:,1}] = deal([1, 2]);
                    [dirs{nRes}]  = deal(vSimp(nRes));
                    [dirs{end}]   = deal(vSimp(end));
        
                    xMax    = max(obj.mesh.coord(:,1));
                    yMax    = max(obj.mesh.coord(:,2));
                    xMin    = min(obj.mesh.coord(:,1));
                    yMin    = min(obj.mesh.coord(:,2));
                    tol     = 1e-10;
        
                    corner1 = @(coor) abs(coor(:,1)-xMin) <= tol & abs(coor(:,2)-yMin)<= tol;
                    corner2 = @(coor) abs(coor(:,1)-xMax) <= tol & abs(coor(:,2)-yMin)<= tol;
                    corner3 = @(coor) abs(coor(:,1)-xMin) <= tol & abs(coor(:,2)-yMax)<= tol;
                    corner4 = @(coor) abs(coor(:,1)-xMax) <= tol & abs(coor(:,2)-yMax)<= tol;
        
                    cornerVec = {corner1; corner2; corner3; corner4};
        
                    isForce = @(coor) (abs(coor(:,1) - xMin)   < 1e-10);
        
        
                    sDir{1}.domain    = @(coor) cornerVec{1}(coor);
                    sDir{1}.direction = cell2mat(dirs(1));
                    sDir{1}.value     = 0;
        
                
                    sDir{2}.domain    = @(coor) cornerVec{2}(coor);
                    sDir{2}.direction = cell2mat(dirs(2));
                    sDir{2}.value     = 0;
        
                    sDir{3}.domain    = @(coor) cornerVec{3}(coor);
                    sDir{3}.direction = cell2mat(dirs(3));
                    sDir{3}.value     = 0;
        
                    sDir{4}.domain    = @(coor) cornerVec{4}(coor);
                    sDir{4}.direction = cell2mat(dirs(4));
                    sDir{4}.value     = 0;
        
                    sDir{5}.domain    = @(coor) cornerVec{nRes}(coor);
                    sDir{5}.direction = cell2mat(dirs(end));
                    sDir{5}.value     = 1;
        
                    % sPL{1}.domain    = @(coor) isForce(coor);
                    % sPL{1}.direction = 2;
                    % sPL{1}.value     = 0;
                    sPL = {};
        
                    dirichletFun = [];
                    for i = 1:numel(sDir)
                        dir = DirichletCondition(obj.mesh, sDir{i});
                        dirichletFun = cat(2, dirichletFun, dir);
                    end
                    s.dirichletFun = dirichletFun;
        
                    pointloadFun = [];
                    for i = 1:numel(sPL)
                        pl = PointLoad(obj.mesh, sPL{i});
                        pointloadFun = cat(2, pointloadFun, pl);
                    end
                    s.pointloadFun = pointloadFun;
        
                    s.periodicFun  = [];
                    s.mesh = obj.mesh;
                    bc = BoundaryConditions(s);
            
                end

            function createDisplacementFunHere(obj) % Crea Funcio FEM pels elements finits
                obj.displacementFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            end

            function createBCApplyerHere(obj, cParams)
                s.mesh = obj.mesh;
                s.boundaryConditions = cParams.boundaryConditions;
                obj.bcApplier = BCApplier(s);
            end

            function createSolverHere(obj,cParams)  % Especifica els parametres del solver
                sS.type      = cParams.solverCase;
                solver       = Solver.create(sS);
                p.solverType = cParams.solverType; % e.g 'Reduced' 'monolithic'
                p.solverMode = cParams.solverMode; %'DIRECT'
                p.solver     = solver;    
%                p.boundaryConditions = cParams.boundaryConditions;
                p.BCApplier          = obj.bcApplier;
%                obj.problemSolver    = ProblemSolver(p);
            end


            function computeStiffnessMatrixHere(obj) % Per obtenir matriu K   
                C     = obj.material;
                f = @(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u)));
                obj.stiffness = IntegrateLHS(f,obj.displacementFun,obj.displacementFun,obj.mesh,'Domain',2);
            end

            function computeForcesHere(obj,cParams)   % Per obtenir vector F
                n.type     = 'Elastic';
                n.scale    = 'MACRO';
                n.dim      = obj.getFunDimsHere();
               % n.BC       = cParams.boundaryConditions;
                n.mesh     = obj.mesh;
                n.material = obj.material;
                n.globalConnec = obj.mesh.connec;
                RHSint = RHSIntegrator.create(n);
                rhs = RHSint.compute();
                % Perhaps move it inside RHSint?
                if strcmp(cParams.solverType,'REDUCED')
                    R = RHSint.computeReactions(obj.stiffness);
                    obj.forces = rhs+R;
                else
                    obj.forces = rhs;
                end

            end

                function dim = getFunDimsHere(obj)
                    d.ndimf     = obj.displacementFun.ndimf;
                    d.nnodes    = size(obj.displacementFun.fValues, 1);
                    d.ndofs     = d.nnodes*d.ndimf;
                    d.nnodeElem = obj.mesh.nnodeElem; % should come from interp..
                    d.ndofsElem = d.nnodeElem*d.ndimf;
                    dim         = d;
                end

            function Cg = computeCmatP1(obj)  % Per obtenir la c
                s.quadType = 2;
                s.boundaryMeshJoined    = obj.boundaryMeshJoined; %Uneix els 4 edges de la mesh
                s.localGlobalConnecBd   = obj.localGlobalConnecBd;
                s.nnodes                 = obj.mesh.nnodes;
    
                % lhs = LHSintegrator_ShapeFunction_fun(s);
                lhs = LHSintegrator_MassBoundary_albert(s);
                test   = LagrangianFunction.create(obj.boundaryMeshJoined, obj.mesh.ndim, 'P1'); % !!
                 obj.dLambda  = LagrangianFunction.create(obj.boundaryMeshJoined, obj.mesh.ndim, 'P1');
                 Cg = lhs.compute(obj.dLambda,test);      
            end


            function rDir = RHSdir(obj)  % Per obtenir la ud
                test   = LagrangianFunction.create(obj.boundaryMeshJoined, obj.mesh.ndim, 'P1');
                s.mesh = obj.boundaryMeshJoined;
                s.quadType = 2;
                rhs = RHSIntegratorShapeFunction(s);
            
                f1x = @(x) [1/(4)*(1-x(1,:,:)).*(1-x(2,:,:));...
                        0*x(2,:,:)  ];
                f2x = @(x) [1/(4)*(1+x(1,:,:)).*(1-x(2,:,:));...
                        0*x(2,:,:)  ];
                f3x = @(x) [1/(4)*(1+x(1,:,:)).*(1+x(2,:,:));...
                        0*x(2,:,:)  ];
                f4x = @(x) [1/(4)*(1-x(1,:,:)).*(1+x(2,:,:));...
                        0*x(2,:,:)  ];
            
                f1y = @(x) [0*x(1,:,:);...
                        1/(4)*(1-x(1,:,:)).*(1-x(2,:,:))  ];
                f2y = @(x) [0*x(1,:,:);...
                        1/(4)*(1+x(1,:,:)).*(1-x(2,:,:))  ];
                f3y = @(x) [0*x(1,:,:);...
                        1/(4)*(1+x(1,:,:)).*(1+x(2,:,:))  ];
                f4y = @(x) [0*x(1,:,:);...
                        1/(4)*(1-x(1,:,:)).*(1+x(2,:,:))  ];
            
            
                f     = {f1x f1y f2x f2y f3x f3y f4x f4y}; %
                nfun = size(f,2);
                rDir = [];
                for i=1:nfun
                    Ud{i}  = AnalyticalFunction.create(f{i},obj.mesh.ndim,obj.boundaryMeshJoined)
                
                    % %% Project to P1
                    % obj.dLambda{i} = obj.dLambda{i}.project('P1');
            
                    rDire = rhs.compute(Ud{i},test);
                      %[iLoc,jLoc,vals] = find(Ce);
            
                      %  l2g_dof = ((obj.localGlobalConnecBd*test.ndimf)' - ((test.ndimf-1):-1:
                      %  l2g_dof = l2g_dof(:);
                      %  jGlob = l2g_dof(jLoc);
                      %  Cg = [Cg sparse(iLoc,jGlob,vals, obj.displacementFun.nDofs, dLambda.nD
            
                      %l2g_dof = ((obj.localGlobalConnecBd*test.ndimf)' - ((test.ndimf-1):-1:0)
                      %l2g_dof = l2g_dof(:);
                      %iGlob = l2g_dof(iLoc);
                    rDir = [rDir rDire];
                end
            end

            function [u, L] = computeDisplacementHere(obj,c, rdir)
                K = obj.stiffness;
                nC  = size(c,2);
                Z   = zeros(nC);
                LHS = [K, c; c' Z];
                % ud = eye(nC);
% %                 ud(7) = 1;
%                 ud(1) = 1;
                forces = zeros(obj.displacementFun.nDofs, size(rdir, 2));
    
                RHS = [forces; rdir];
                sol = LHS\RHS;
                u = sol(1:obj.displacementFun.nDofs,:);
                L = -sol(obj.displacementFun.nDofs+1:end,:); 
%                 obj.displacementFun.fValues = u;
%                 EIFEMtesting.plotSolution(u,obj.mesh,1,1,0,0)
            end


    end
end
