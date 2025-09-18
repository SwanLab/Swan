classdef ElasticProblemMicro < handle
    
    properties (Access = public)
        uFluc, strain, stress
        Chomog
    end

    properties (Access = private)
        mesh
        material
        trialFun
        boundaryConditions, bcApplier
        solverType, solverMode, solverCase
        lagrangeMultipliers
        problemSolver
    end

    methods (Access = public)

        function obj = ElasticProblemMicro(cParams)
            obj.init(cParams);
            obj.createTrialFun();
            obj.createBCApplier();
            obj.createSolver();
        end

        function obj = solve(obj)
            C     = obj.material;
            f     = @(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u)));
            LHS   = IntegrateLHS(f,obj.trialFun,obj.trialFun,obj.mesh,2);
            for iB = 1:obj.computeNbasis()
                [eB,v] = obj.createDeformationBasis(iB);
                f = @(v) -DDP(SymGrad(v),DDP(C,eB));
                RHS = IntegrateRHS(f,obj.trialFun,obj.mesh,2);    
                uF{iB}      = obj.computeDisplacement(LHS,RHS,iB);
                strainF{iB} = eB+SymGrad(uF{iB});
                stressF{iB} = DDP(obj.material, strainF{iB});
                ChiB        = Integrator.compute(stressF{iB},obj.mesh,2);
                obj.convertChomogToFourthOrder(ChiB,v,iB);
            end
            obj.uFluc  = uF;
            obj.strain = strainF;
            obj.stress = stressF;
        end

        function v = computeGeometricalVolume(obj)
            v = 1;%sum(sum(obj.geometry.dvolu));
        end

        function setMatProps(obj,s)
           obj.material.compute(s);
        end
        
        function updateMaterial(obj, mat)
            obj.material = mat;
        end

        function [fun, funNames] = getFunsToPlot(obj)
            fun = {obj.uFluc, obj.strain.project('P1'), ...
                obj.stress.project('P1')};
            funNames = {'displacement', 'strain', 'stress'};
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh     = cParams.mesh;
            obj.material = cParams.material;
            obj.solverType = cParams.solverType;
            obj.solverMode = cParams.solverMode;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.solverCase  = cParams.solverCase;
        end

        function createTrialFun(obj)
            obj.trialFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
        end

       function [s,v] = createDeformationBasis(obj,iBasis)
           v      = obj.computeBasesPosition();
           sV     = zeros(obj.mesh.ndim,obj.mesh.ndim);
           sV(v(iBasis,1),v(iBasis,2)) = 1;
           sHV = diag(diag(sV));
           sDV = sV-sHV;
           sV = sHV+1*(sDV+sDV');
           s = ConstantFunction.create(sV,obj.mesh);
       end

       function v = computeBasesPosition(obj)
           switch obj.mesh.ndim
               case 2
                   v = [1,1; 2,2; 1,2];
               case 3
                   v = [1,1; 2,2; 3,3; 2,3; 1,3; 1,2];
           end
       end

        function nBasis = computeNbasis(obj)
            homogOrder = 1;
            nDim = obj.mesh.ndim;
            nBasis = homogOrder*nDim*(nDim+1)/2;
        end

        function dim = getFunDims(obj)
            d.ndimf  = obj.trialFun.ndimf;
            d.nnodes = size(obj.trialFun.fValues, 1);
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = obj.mesh.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

        function createBCApplier(obj)
            s.mesh = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            bc = BCApplier(s);
            obj.bcApplier = bc;
        end

        function createSolver(obj)
            sS.type =  obj.solverCase;
            solver = Solver.create(sS);
            s.solverType = obj.solverType;
            s.solverMode = obj.solverMode;
            s.solver     = solver;
            s.boundaryConditions = obj.boundaryConditions;
            s.BCApplier = obj.bcApplier;
            obj.problemSolver = ProblemSolver(s);
        end

        function uFun = computeDisplacement(obj, LHS, RHS, iB)
            s.stiffness = LHS;
            s.forces    = RHS;
            s.iBase     = iB;
            s.nBasis    = obj.computeNbasis();
            [u, L]      = obj.problemSolver.solve(s);
            obj.lagrangeMultipliers = L;
            uSplit = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            uFun = copy(obj.trialFun);
            uFun.setFValues(uSplit);
        end

        function convertChomogToFourthOrder(obj,ChiB,v,iB)
            Ch = obj.Chomog;
            v1 = v(iB,1);    v2 = v(iB,2);
            if v1==v2
                Ch(:,:,v1,v2) = ChiB;
            else
                ChShear        = zeros(size(ChiB));
                ChShear(v1,v2) = ChiB(v1,v2);
                Ch(:,:,v1,v2)  = ChShear;
                ChShear        = zeros(size(ChiB));
                ChShear(v2,v1) = ChiB(v2,v1);
                Ch(:,:,v2,v1)  = ChShear;
            end
            obj.Chomog = Ch;
        end

    end

    

end
