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
            ndim = obj.mesh.ndim;
            LHS = obj.computeLHS();
            nBasis = obj.computeNbasis();
            obj.Chomog = zeros(ndim, ndim, ndim, ndim);
            for iB = 1:nBasis
                [strainB,v] = obj.createDeformationBasis(iB);
                RHS         = obj.computeRHS(strainB,LHS);
                uF{iB}      = obj.computeDisplacement(LHS,RHS,iB,nBasis);
                strainF{iB} = strainB+SymGrad(uF{iB});
                stressF{iB} = DDP(obj.material, strainF{iB});
                ChiB        = obj.computeChomog(stressF{iB});
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

        function LHS = computeLHS(obj)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.quadratureOrder = 2;
            lhs = LHSIntegrator.create(s);
            C     = obj.material;
            for i = 1:obj.trialFun.nDofsElem
                v = Test(obj.trialFun,i);
                for j = 1:obj.trialFun.nDofsElem
                    u = Test(obj.trialFun,j);
                    f{i,j} = DDP(SymGrad(v),DDP(C,SymGrad(u)));
                end
            end
            LHS = lhs.compute(f,obj.trialFun,obj.trialFun);

        end

        function rhs = computeRHS(obj,strainBase,LHS)
            s.fun  = obj.trialFun;
            s.type = 'ElasticMicro';
            s.dim      = obj.getFunDims();
            s.BC       = obj.boundaryConditions;
            s.mesh     = obj.mesh;
            s.material = obj.material;
            s.globalConnec = obj.mesh.connec;
            RHSint = RHSIntegrator.create(s);

            C = obj.material;
            for i = 1:obj.trialFun.nDofsElem
                v = Test(obj.trialFun,i);
                f{i} = DDP(SymGrad(v),DDP(C,strainBase));
            end
            rhs = RHSint.compute(f,obj.trialFun);
            R = RHSint.computeReactions(LHS); %%?
        end

        function uFun = computeDisplacement(obj, LHS, RHS, iB, nBasis)
            s.stiffness = LHS;
            s.forces    = RHS;
            s.iBase     = iB;
            s.nBasis    = nBasis;
            [u, L]      = obj.problemSolver.solve(s);
            obj.lagrangeMultipliers = L;
            uSplit = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            uFun = copy(obj.trialFun);
            uFun.setFValues(uSplit);
        end

        function Chomog = computeChomog(obj,stress)
            Chomog = Integrator.compute(stress,obj.mesh,2);
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
