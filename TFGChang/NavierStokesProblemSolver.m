classdef NavierStokesProblemSolver < handle

    properties (Access = public)
        velocityFun
        pressureFun
    end

    properties (Access = private)
        mesh
        forcesFormula
        dirConditions
        dirDofs
        material
        velocityField
        nu
        dt
    end

    properties (Access = private)
        diffTime
        numDofs
        solver
        LHSS
        LHS
        RHS
        x
        fullX
        u
        p
        M
        D
        KU
        KP
        C
        uT
    end

    methods (Access = public)

        function obj = NavierStokesProblemSolver(cParams)
            obj.init(cParams);
            obj.setDiffTime();
            obj.calculateNumDofs();
            obj.createSolver();
        end

        function compute(obj)
            obj.solveNVProblem();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh             = cParams.mesh;
            obj.velocityFun      = cParams.velocityFun;
            obj.pressureFun      = cParams.pressureFun;
            obj.forcesFormula    = cParams.forcesFormula;
            obj.dirConditions    = cParams.dirConditions;
            obj.dirDofs          = cParams.dirDofs;
            obj.material         = cParams.material;
            obj.velocityField    = cParams.velocityField;
            obj.nu               = cParams.nu;
            obj.dt               = cParams.dt;
        end

        function setDiffTime(obj)
            obj.diffTime = Inf;
        end

        function calculateNumDofs(obj)
            obj.numDofs = obj.velocityFun.nDofs + obj.pressureFun.nDofs;
        end

        function createSolver(obj)
            b.type =  'DIRECT';
            obj.solver = Solver.create(b);
        end

        function computeLHSStokes(obj)
            c.type          = 'Stokes';
            c.dt            = obj.diffTime;
            c.mesh          = obj.mesh;
            c.material      = obj.material;
            c.velocityFun   = obj.velocityFun;
            c.pressureFun   = obj.pressureFun;
            LHSintegratorg  = LHSIntegrator.create(c);
            obj.LHS         = LHSintegratorg.compute();
            obj.LHSS        = obj.LHS;
        end

        function computeLHSNS(obj,residual)
            c.type          = 'NavierStokes';
            c.mesh          = obj.mesh;
            c.material      = obj.material;
            c.velocityField = obj.velocityField;
            c.velocityFun   = obj.velocityFun;
            c.pressureFun   = obj.pressureFun;
            c.LHSS          = obj.LHSS;
            c.dt            = obj.dt;
            c.residual      = residual;
            LHSintegratorg  = LHSIntegrator.create(c);
            obj.LHS         = LHSintegratorg.compute();
        end

        function computeLHSMatrixNS(obj)
            c.type                          = 'NavierStokes';
            c.mesh                          = obj.mesh;
            c.material                      = obj.material;
            c.velocityField                 = obj.velocityField;
            c.velocityFun                   = obj.velocityFun;
            c.pressureFun                   = obj.pressureFun;
            c.dt                            = obj.dt;
            LHSintegratorg                  = LHSIntegrator.create(c);
            [obj.M, obj.D, obj.KU, obj.KP]  = LHSintegratorg.computeLHSMatrix();
        end

        function computeLHSConvectiveMatrix(obj)
            c.type          = 'NavierStokes';
            c.mesh          = obj.mesh;
            c.material      = obj.material;
            c.velocityField = obj.velocityField;
            c.velocityFun   = obj.velocityFun;
            c.dt            = obj.dt;
            LHSintegratorg  = LHSIntegrator.create(c);
            obj.C           = LHSintegratorg.computeConvectiveTerm();
        end

        function computeTentativeVelocity(obj)
            uN  = obj.velocityField.fValues;
            uN  = reshape(uN.', [], 1);
            LHSTV   = obj.M + obj.C + obj.KU;
            RHSTV   = obj.M * uN;
            obj.uT  = obj.solver.solve(LHSTV, RHSTV);
        end

        function computePressure(obj)
            RHSP     = obj.D * obj.uT ./ obj.dt;
            obj.p    = obj.solver.solve(obj.KP,RHSP);
        end

        function computeVelocity(obj)
            obj.u = obj.uT - obj.dt * obj.D * obj.p;
        end

        function computeRHS(obj)
            d.type          = 'Stokes';
            d.mesh          = obj.mesh;
            d.velocityFun   = obj.velocityFun;
            d.pressureFun   = obj.pressureFun;
            d.forcesFormula = obj.forcesFormula;
            RHSint          = RHSIntegrator.create(d);
            F               = RHSint.integrate();
            uD      = obj.dirConditions(:,3);
            R       = -obj.LHS(:,obj.dirDofs)*uD;
            obj.RHS = F + R;
        end

        function computeSolution(obj)
            freeDofsPlus = setdiff(1:obj.numDofs,obj.dirDofs);
            LHSr         = obj.LHS(freeDofsPlus,freeDofsPlus);
            RHSr         = obj.RHS(freeDofsPlus);
            obj.x        = obj.solver.solve(LHSr, RHSr);
        end

        function addDirBoundaryConditions(obj)
            uD                 = obj.dirConditions(:,3);
            nSteps             = length(obj.x(1,:));
            uD                 = repmat(uD,1,nSteps);
            obj.fullX          = zeros(obj.numDofs,nSteps);
            freeDofs           = setdiff(1:(obj.numDofs),obj.dirDofs);
            obj.fullX(freeDofs,:) = obj.x;
            if ~isempty(obj.dirDofs)
                obj.fullX(obj.dirDofs,:) = uD;
            end
        end

        function addDirBCToVelocity(obj)
            uD                 = obj.dirConditions(1:end - 1,3);
            uDirDofs           = obj.dirDofs(1:end - 1);
            if ~isempty(uDirDofs)
                obj.u(uDirDofs,:) = uD;
            end
        end

        function addDirBCToPressure(obj)
            pD                 = obj.dirConditions(end,3);
            pDirDofs           = obj.dirDofs(end);
            if ~isempty(pDirDofs)
                obj.p(pDirDofs,:) = pD;
            end
        end

        function separateVariables(obj)
            ndofsV = obj.velocityFun.nDofs;
            obj.u = obj.fullX(1:ndofsV,:);
            obj.p = obj.fullX(ndofsV+1:end,:);
        end

        function defineVariablesFValues(obj)
            nFieldu = obj.velocityFun.ndimf;
            nnode   = round(length(obj.u)/nFieldu);
            nodes   = 1:nnode;
            velfval = zeros(nnode,nFieldu);
            for idim = 1:nFieldu
                dofs = nFieldu*(nodes-1)+idim;
                velfval(:,idim) = obj.u(dofs, end);
            end
            obj.velocityFun.setFValues(velfval);
            obj.pressureFun.setFValues(obj.p(:,end));
        end

        function initPicardIteration(obj)
            obj.computeLHSStokes();
            %obj.computeLHSMatrixNS();
            %obj.computeLHSNS();
            obj.computeRHS();
            obj.computeSolution();
            obj.addDirBoundaryConditions();
            obj.separateVariables();
            obj.defineVariablesFValues();
        end

        function initPicardIterationFSM(obj)
            obj.computeLHSStokes();
            obj.computeLHSMatrixNS();
            %obj.computeLHSNS();
            obj.computeRHS();
            obj.computeSolution();
            obj.addDirBoundaryConditions();
            obj.separateVariables();
            obj.defineVariablesFValues();
        end

        function computePicardIteration(obj,residual)
            obj.computeLHSNS(residual);
            %obj.computeLHSMatrixNS();
            obj.computeRHS();
            obj.computeSolution();
            obj.addDirBoundaryConditions();
            obj.separateVariables();
            obj.defineVariablesFValues();
        end

        function [Residual, ReError] = computeError(obj)
            vNew = obj.velocityFun.fValues;
            vOld = obj.velocityField.fValues;
            denom = max(abs(vNew), 1e-12);

            Residual = max(max(abs(vNew - vOld)));
            ReError  = max(max(abs((vNew - vOld) ./ denom)));
        end

        function solveNVProblemFSM(obj)

            Residual     = 1;
            ReError      = 1;
            iter         = 1;
            maxIter      = 1e3;
            tolRes       = 1e-6;
            tolReE       = 1e-4;

            if obj.nu < 1e-3
                relaxFactor = 0.5;
            elseif obj.nu < 1e-2
                relaxFactor = 0.7;
            else
                relaxFactor = 1;
            end

            obj.initPicardIterationFSM()

            fprintf('   Picard Iteration             Residual                ReError\n');

            while ~((Residual <= tolRes && ReError <= tolReE) || iter >= maxIter)

                obj.computeLHSConvectiveMatrix();
                obj.computeTentativeVelocity();
                obj.computePressure();
                obj.computeVelocity();

                obj.addDirBCToVelocity();
                obj.defineVariablesFValues();

                [Residual, ReError] = obj.computeError();

                fprintf(['      %d               ' ...
                    '        %.7e          %.5e\n'], iter, Residual, ReError);

                relaxedVel = (1 - relaxFactor) * obj.velocityField.fValues + relaxFactor * obj.velocityFun.fValues;
                obj.velocityField.setFValues(relaxedVel);

                iter = iter + 1;

            end

        end

        function solveNVProblem(obj)

            Residual     = 1;
            ReError      = 1;
            iter         = 1;
            maxIter      = 1e3;
            tolRes       = 1e-6;
            tolReE       = 1e-4;

            if obj.nu < 1e-3
                relaxFactor = 0.5;
            elseif obj.nu < 1e-2
                relaxFactor = 0.7;
            else
                relaxFactor = 1;
            end

            if all(obj.velocityField.fValues(:) == 0)
                obj.initPicardIteration();
                obj.velocityField.setFValues(obj.velocityFun.fValues);
            else

                obj.computeLHSStokes()
            end
            fprintf('   Picard Iteration             Residual                ReError\n');

            while ~((Residual <= tolRes && ReError <= tolReE) || iter >= maxIter)

                obj.computePicardIteration(Residual);

                [Residual, ReError] = obj.computeError();

                fprintf(['      %d               ' ...
                    '        %.7e          %.5e\n'], iter, Residual, ReError);

                relaxedVel = (1 - relaxFactor) * obj.velocityField.fValues + relaxFactor * obj.velocityFun.fValues;
                obj.velocityField.setFValues(relaxedVel);

                iter = iter + 1;

            end

        end


    end

end