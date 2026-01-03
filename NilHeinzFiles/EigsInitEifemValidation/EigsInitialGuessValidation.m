classdef EigsInitialGuessValidation < handle

    properties (Access = public)

    end
    properties (Access = private)
        referenceMesh
        meshDomain
        boundaryConditions
        bcApplier
        LHS
        RHS
        fileNameEIFEM
        tolSameNode
        solverType
    end


    properties (Access = private)
        nSubdomains
    end

    methods (Access = public)

        function obj = EigsInitialGuessValidation()
            close all
            obj.init()
            
            %Mesh
            radiusMesh = 0.25;
            obj.createReferenceMesh(false,radiusMesh);
            bS  = obj.referenceMesh.createBoundaryMesh();
            obj.meshDomain = obj.referenceMesh;
            
            %Boundary Conditions
            [bC,dir] = obj.createBoundaryConditions();
            obj.boundaryConditions = bC;
            obj.createBCapplier()

            %K and M generation
            [LHSr, Mr, RHSr] = obj.createElasticProblem();
           
            %% Modal Analysis
            
            


            % Assume you already have the true first eigenvector in X_true
            % and matrices K and M (stiffness and mass, or whatever your generalized problem is)
            % e.g. from a full eig solve once at the beginning (expensive but done only once)
            
            % Example (uncomment if you need to compute the reference once):
            % [X_full, D_full] = eig(K, M, 'chol');
            % [lambda_full, idx] = sort(diag(D_full));
            % X_full = X_full(:, idx);
            % X_true = X_full(:, 1);         % true dominant/smallest eigenvector
            % true_eigval = lambda_full(1);
            
           % 1. SETUP: Get your system matrices and reference solution
            % Assuming LHSr is Stiffness (K) and Mr is Mass (M)
            K = LHSr; 
            M = Mr;
            
            % Get the "Truth" (Reference)
            % We extract just the 1st column to ensure X_true is a vector, not a matrix
            [true_eigenvalues, true_eigenvectors, natFreq] = obj.computeModalAnalysis(LHSr, Mr);
            
                       
            X_true = true_eigenvectors(:, 1); % Select the specific mode you are targeting
            true_eigval = true_eigenvalues(1);
            
            % 2. PARAMETERS
            eigs_to_get = 1;               
            sigma_values = logspace(0, 4, 30); 
            nSigma = length(sigma_values);
            
            % Preallocate
            times     = zeros(nSigma, 1);
            iterations = zeros(nSigma, 1);
            computed_eigval = zeros(nSigma, 1);
            residual_norms  = zeros(nSigma, 1);
            
            % Base options
            % Note: We set a tight tolerance. If sigma is 1e-12 but tol is 1e-6, 
            % eigs will exit immediately for all small sigmas, hiding the behavior you want to study.
            opts = struct('disp', 2, 'tol', 1e-14); 
            
            fprintf('sigma          time (s)      iterations    computed lambda      rel. error\n');
            fprintf('%s\n', repmat('-', 1, 95));
            
            for i = 1:nSigma
                sigma = sigma_values(i);
                
                % 3. PERTURBATION
                rng(12345 + i); 
                noise = randn(size(X_true));
                noise = noise / norm(noise);    % Unit norm noise
                
                % Create guess: Truth + Sigma * Noise
                v0 = X_true + sigma * noise;
                v0 = v0 / norm(v0);             % Normalize guess
                
                % Set StartVector
                opts.v0 = v0; 
                opts.p = 3;
                
                % 4. SOLVE WITH EVALC
                % We use evalc to capture the command window text to parse iteration counts.
                tic;
                % Note: variables V_calc, D_calc, flag must be uniquely named to avoid conflicts
                txt_output = evalc('[V_calc, D_calc, flag] = eigs(K, M, eigs_to_get, "smallestabs", opts);');
                %[Vb,Db] = eigs(K,M,1,'smallestabs',opts);
                elapsed = toc;
                
                % 5. PARSE ITERATIONS
                % Regex to find "Iteration X: ..." patterns standard in ARPACK output
                iter_matches = regexp(txt_output, 'Iteration\s+(\d+)', 'tokens');

                if ~isempty(iter_matches)
                    % Convert cell strings to numbers and find the max
                    iter_nums = cellfun(@str2double, [iter_matches{:}]);
                    iterations(i) = max(iter_nums);
                else
                    % If no text, it likely converged instantly (0 iterations usually means 
                    % the starting vector was within tolerance immediately)
                    iterations(i) = 0; 
                end

                % 6. STORE RESULTS
                times(i) = elapsed;

                % Extract scalar eigenvalue
                if isempty(D_calc)
                    computed_eigval(i) = NaN; % Handle non-convergence
                else
                    computed_eigval(i) = D_calc(1,1);

                    % Calculate Residual: ||KV - M*V*lambda|| / ||KV||
                    % (Using V_calc to ensure dimensions match)
                    res_vec = K*V_calc - M*V_calc*D_calc(1,1);
                    residual_norms(i) = norm(res_vec) / norm(K*V_calc);
                end

                % Display
                rel_err = abs(computed_eigval(i) - true_eigval)/abs(true_eigval);
                fprintf('%6.2e     %8.4f       %6d        %.12e     %.2e\n', ...
                        sigma, elapsed, iterations(i), computed_eigval(i), rel_err);
            end
            
            % 7. PLOT RESULTS (Optional)
            figure;
            semilogx(sigma_values, iterations, '-o');
            xlabel('Sigma (Perturbation)');
            ylabel('Iterations to Converge');
            title('Effect of Initial Guess Quality on Eigs Convergence');
            grid on;
            
            %sgtitle('ARPACK/eigs convergence study (generalized eigenvalue problem)');
        
            
        end

    end

    methods (Access = private)

        function init(obj)
            obj.nSubdomains  = [1 1]; %nx ny
            obj.fileNameEIFEM = 'DEF_Q4porL_1.mat';
            obj.tolSameNode = 1e-10;
            obj.solverType = 'REDUCED';
        end        

        function createReferenceMesh(obj, loadfile, radius)

            holeMesh    = obj.createMesh(radius);
            s.coord     = holeMesh.coord;
            s.connec    = holeMesh.connec;
            s.interType = 'LINEAR';
            s           = obj.updateCoordsMesh(s); 
            obj.referenceMesh = Mesh.create(s);
        end

        function s = updateCoordsMesh(obj, s)
            % Nudge nodes at the four rectangle corners in x to avoid
            % exact coincidences.
            tol  = 1e-8;
            epsx = 1e-9;
        
            x = s.coord(:,1); y = s.coord(:,2);
            xmax = max(x); xmin = min(x);
            ymax = max(y); ymin = min(y);
        
            % Top-right (xmax,ymax)
            mask = abs(x - xmax) < tol & abs(y - ymax) < tol;
            s.coord(mask, :) = s.coord(mask, :) - [epsx, 0];
        
            % Bottom-right (xmax,ymin)
            mask = abs(x - xmax) < tol & abs(y - ymin) < tol;
            s.coord(mask, :) = s.coord(mask, :) - [epsx, 0];
        
            % Top-left (xmin,ymax)
            mask = abs(x - xmin) < tol & abs(y - ymax) < tol;
            s.coord(mask, :) = s.coord(mask, :) + [epsx, 0];
        
            % Bottom-left (xmin,ymin)
            mask = abs(x - xmin) < tol & abs(y - ymin) < tol;
            s.coord(mask, :) = s.coord(mask, :) + [epsx, 0];
        end
        
        function holeMesh = createMesh(obj,radius)
            fullmesh = UnitTriangleMesh(720,720);
            ls = obj.computeCircleLevelSet(fullmesh,radius);
            sUm.backgroundMesh = fullmesh;
            sUm.boundaryMesh   = fullmesh.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(ls);
            holeMesh = uMesh.createInnerMesh();
            
        end
  
        function ls = computeCircleLevelSet(obj, mesh,radius)
            gPar.type          = 'Circle';
            gPar.radius        = radius;
            gPar.xCoorCenter   = 0.5;
            gPar.yCoorCenter   = 0.5;
            g                  = GeometricalFunction(gPar);
            phiFun             = g.computeLevelSetFunction(mesh);
            lsCircle           = phiFun.fValues;
            ls = -lsCircle;
        end
     

        function [mD,mSb,iC,lG,iCR,discMesh] = createMeshDomain(obj)
            s.nsubdomains   = obj.nSubdomains;
            s.meshReference = obj.referenceMesh;
            s.tolSameNode = obj.tolSameNode;
            m = MeshCreatorFromRVE.create(s);
            [mD,mSb,iC,~,lG,iCR,discMesh] = m.create();
        end

        
        function mCoarse = createCoarseMesh(obj,mR)
            s.nsubdomains   = obj.nSubdomains; %nx ny
            s.meshReference = obj.createReferenceCoarseMesh(mR);
            s.tolSameNode   = obj.tolSameNode;
            mRVECoarse      = MeshCreatorFromRVE.create(s);
            [mCoarse,~,~] = mRVECoarse.create();
        end

        function cMesh = createReferenceCoarseMesh(obj,mR)
            xmax = max(mR.coord(:,1));
            xmin = min(mR.coord(:,1));
            ymax = max(mR.coord(:,2));
            ymin = min(mR.coord(:,2));
            coord(1,1) = xmin;  coord(1,2) = ymin;
            coord(2,1) = xmax;  coord(2,2) = ymin;
            coord(3,1) = xmax;  coord(3,2) = ymax;
            coord(4,1) = xmin;  coord(4,2) = ymax;
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
            E  = 1;
            nu = 1/3;  
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = mesh.ndim;
            s.young   = ConstantFunction.create(E,mesh);
            s.poisson = ConstantFunction.create(nu,mesh);
            tensor    = Material.create(s);
            material  = tensor;
        end


        function [bC,Dir] = createBoundaryConditions(obj)
            minx = min(obj.meshDomain.coord(:,1));
            maxx = max(obj.meshDomain.coord(:,1));
            tolBound = obj.tolSameNode;
            isLeft   = @(coor) (abs(coor(:,1) - minx)   < tolBound);
            isRight  = @(coor) (abs(coor(:,1) - maxx)   < tolBound);
            Dir{1}.domain    = @(coor) isLeft(coor);%| isRight(coor) ;
            Dir{1}.direction = [1,2];
            Dir{1}.value     = 0;
            dirichletFun = DirichletCondition(obj.meshDomain, Dir{1});

            mesh = obj.meshDomain;
            PL.domain    = @(coor) isRight(coor);
            PL.direction = 2;
            PL.value     = -0.1;
            pointload = PointLoad(mesh,PL);
            % need this because force applied in the face not in a point
            pointload.values        = pointload.values/size(pointload.dofs,1);
            fvalues                 = zeros(mesh.nnodes*mesh.ndim,1);
            fvalues(pointload.dofs) = pointload.values;
            fvalues                 = reshape(fvalues,mesh.ndim,[])';
            pointload.fun.setFValues(fvalues);

            s.pointloadFun = pointload;
            s.dirichletFun = dirichletFun;
            s.periodicFun  =[];
            s.mesh         = mesh;
            bC             = BoundaryConditions(s);                        
        end

  

        function [LHSr,Mr, RHSr] = createElasticProblem(obj)
            u = LagrangianFunction.create(obj.meshDomain,obj.meshDomain.ndim,'P1');
            material = obj.createMaterial(obj.meshDomain);
            [lhs,LHSr] = obj.computeStiffnessMatrix(obj.meshDomain,u,material);
            [M,Mr] = obj.computeMassMatrix(obj.meshDomain,u);
            RHSr       = obj.computeForces(lhs,u);
        end

        function [LHS,LHSr] = computeStiffnessMatrix(obj,mesh,dispFun,C)

            LHS = IntegrateLHS(@(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u))),dispFun,dispFun,mesh,2);
            LHSr = obj.bcApplier.fullToReducedMatrixDirichlet(LHS);
        end

        function [M , Mr] = computeMassMatrix(obj, mesh, dispFun)
            rho = obj.computeDensity(mesh);
            M = IntegrateLHS(@(u,v) rho .* DP(v,u),dispFun,dispFun,mesh,2); %DP(u,v) 
            Mr = obj.bcApplier.fullToReducedMatrixDirichlet(M);

        end

        function rho = computeDensity(obj, mesh)
            rho = 1; % kg/m^3 -
            rho  = ConstantFunction.create(rho, mesh);
        end

        function [lambda, Phi, omega] = computeModalAnalysis(obj, K, M)
            [Phi, D] = eigs(K, M, 15, "smallestabs");
            lambda = diag(D);
            [lambda, idx] = sort(lambda, 'ascend'); 
            Phi = Phi(:, idx); %sort
            omega = sqrt(max(lambda,0));
        end

        function RHS = computeForces(obj,stiffness,u)
           
            ndofs = u.nDofs;
            bc            = obj.boundaryConditions;
            neumann       = bc.pointload_dofs;
            neumannValues = bc.pointload_vals;
            rhs = zeros(ndofs,1);
            if ~isempty(neumann)
                rhs(neumann) = neumannValues;
            end
            if strcmp(obj.solverType,'REDUCED')
                bc      = obj.boundaryConditions;
                dirich  = bc.dirichlet_dofs;
                dirichV = bc.dirichlet_vals;
                if ~isempty(dirich)
                    R = -stiffness(:,dirich)*dirichV;
                else
                    R = zeros(sum(ndofs(:)),1);
                end
                rhs = rhs+R;
            end
            RHS = obj.bcApplier.fullToReducedVectorDirichlet(rhs);
        end

        function [Kcoarse, Mcoarse, T] = offlineTraining(obj) %(obj,dir,iC,lG,bS,iCR,dMesh,radiusMesh)
            mR = obj.referenceMesh;
          
            data = Training(mR);
            p = OfflineDataProcessor(data); % i don't want to have to run this if I use NN
            EIFEoper = p.computeROMbasis();
            
            Kcoarse = EIFEoper.Kcoarse;
            Mcoarse = EIFEoper.Mcoarse;    
            Udef = EIFEoper.Udef;
            Urb = EIFEoper.Urb;
            T = Udef + Urb;


            
        end
        
        function d = createDomainDecompositionDofManager(obj,iC,lG,bS,mR,iCR)
            s.nSubdomains     = obj.nSubdomains;
            s.interfaceConnec = iC;
            s.interfaceConnecReshaped = iCR;
            s.locGlobConnec   = lG;
            s.nBoundaryNodes  = bS{1}.mesh.nnodes;
            s.nReferenceNodes = mR.nnodes;
            s.nNodes          = obj.meshDomain.nnodes;
            s.nDimf           = obj.meshDomain.ndim;
            d = DomainDecompositionDofManager(s);
        end

        function Milu = createILUpreconditioner(obj,LHS)
            s.LHS = LHS;
            s.type = 'ILU';
            M = Preconditioner.create(s);
            Milu = @(r) M.apply(r);
        end

        
       
        
    end

end
