classdef StokesProblemSolver < handle
    
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
    end

    properties (Access = private)
        diffTime
        numDofs
        solver
        LHS
        RHS
        x
        fullX
        u
        p
    end

    methods (Access = public)
        
        function obj = StokesProblemSolver(cParams)
            obj.init(cParams);
            obj.setDiffTime();
            obj.calculateNumDofs();
            obj.createSolver();
        end

        function compute(obj)
            obj.computeLHS();
            obj.computeRHS();
            obj.solveProblem();
            obj.addDirBoundaryConditions();
            obj.separateVariablesFvalues();
            obj.defineVariables();
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

        function computeLHS(obj)
            c.type          = 'Stokes';
            c.dt            = obj.diffTime;
            c.mesh          = obj.mesh;
            c.material      = obj.material;
            c.velocityFun   = obj.velocityFun;
            c.pressureFun   = obj.pressureFun;
            LHSintegratorg  = LHSintegrator.create(c);
            obj.LHS         = LHSintegratorg.compute();
        end

        function computeRHS(obj)
            d.type          = 'Stokes';
            d.mesh          = obj.mesh;
            d.velocityFun   = obj.velocityFun;
            d.pressureFun   = obj.pressureFun;
            d.forcesFormula = obj.forcesFormula;
            RHSint          = RHSintegrator.create(d);
            F               = RHSint.integrate();
            uD      = obj.dirConditions(:,3);
            R       = -obj.LHS(:,obj.dirDofs)*uD;
            obj.RHS = F + R;
        end

        function solveProblem(obj)
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

        function separateVariablesFvalues(obj)
            ndofsV = obj.velocityFun.nDofs;
            obj.u = obj.fullX(1:ndofsV,:);
            obj.p = obj.fullX(ndofsV+1:end,:);
        end

        function defineVariables(obj)
            nFieldu = obj.velocityFun.ndimf;
            nnode   = round(length(obj.u)/nFieldu);
            nodes   = 1:nnode;
            velfval = zeros(nnode,nFieldu);
            for idim = 1:nFieldu
                dofs = nFieldu*(nodes-1)+idim;
                velfval(:,idim) = obj.u(dofs, end);
            end
            obj.velocityFun.fValues = velfval;
            obj.pressureFun.fValues = obj.p(:,end);
        end


    end

end