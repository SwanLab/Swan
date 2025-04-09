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
            LHS = obj.computeLHS();
            %    oX     = zeros(obj.getDimensions().ndimf,1);
            nBasis = obj.computeNbasis();
            obj.Chomog = zeros(nBasis, nBasis);
            for iB = 1:nBasis
                strainBase = obj.createDeformationBasis(iB);
                RHS      = obj.computeRHS(strainBase,LHS);
                uF{iB}    = obj.computeDisplacement(LHS,RHS);
                strnFluc = SymGrad(uF{iB});
                str{iB} = obj.computeStrain(strainBase,strnFluc);
                sigma{iB} = DDP(obj.material, str{iB});
                sigma{iB}.ndimf = str{iB}.ndimf;
                Ch(:,iB) = obj.computeChomog(sigma{iB});
            end
            obj.uFluc = uF;
            obj.strain = str;
            obj.stress = sigma;
            obj.Chomog = Ch;
        end

        function s = createDeformationBasis(obj,iBasis)
            nBasis = obj.computeNbasis();
            sV = zeros(nBasis,1);
            sV(iBasis) = 1;
            s = ConstantFunction.create(sV,obj.mesh);
            %fH = @(xV) x(1,:,:);
            %x1 = AnalyticalFunction.create(fH,obj.mesh)
        end

        function nBasis = computeNbasis(obj)
            homogOrder = 1;
            nDim = obj.mesh.ndim;
            nBasis = homogOrder*nDim*(nDim+1)/2;
        end

        function LHS = computeLHS(obj)
            ndimf = obj.trialFun.ndimf;
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.test     = LagrangianFunction.create(obj.mesh,ndimf, 'P1');
            s.trial    = obj.trialFun;
            s.material = obj.material;
            s.quadratureOrder = 2;
            lhs = LHSIntegrator.create(s);
            LHS = lhs.compute();
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

        function rhs = computeRHS(obj,strainBase,LHS)
            s.fun  = obj.trialFun;
            s.type = 'ElasticMicro';
            s.dim      = obj.getFunDims();
            s.BC       = obj.boundaryConditions;
            s.mesh     = obj.mesh;
            s.material = obj.material;
            s.globalConnec = obj.mesh.connec;
            RHSint = RHSIntegrator.create(s);
            rhs = RHSint.compute(strainBase);
            R = RHSint.computeReactions(LHS); %%?
        end

        function uFun = computeDisplacement(obj, LHS, RHS)
            s.stiffness = LHS;
            s.forces    = RHS;
        %    s.iVoigt    = iVoigt;
        %    s.nVoigt    = size(obj.forces,2);
            [u, L]      = obj.problemSolver.solve(s);
            obj.lagrangeMultipliers = L;

            uSplit = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            uFun = copy(obj.trialFun);
            uFun.setFValues(uSplit);            

            % z.mesh    = obj.mesh;
            % z.fValues = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            % z.order   = 'P1';
            % uFeFun = LagrangianFunction(z);
            % obj.uFun{iVoigt} = uFeFun;

          %  uSplit = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
          %  obj.displacementFun.setFValues(uSplit);
        end

        function strainFun = computeStrain(obj, strainBase,strainFluc)
            s.operation    = @(xV) strainBase.evaluate(xV)+strainFluc.evaluate(xV);
            s.mesh         = obj.mesh;
            s.ndimf        = size(strainBase,1);
            strainFun      = DomainFunction(s);
        end

        %% 

        function Chomog = computeChomog(obj,stress)
            Chomog = Integrator.compute(stress,obj.mesh,2);
        end

        function Chomog = computeChomogFromLagrangeMultipliers(obj,iBase)
           if strcmp(obj.solverMode, 'DISP')
                L = obj.lagrangeMultipliers;
                nPeriodic = length(obj.boundaryConditions.periodic_leader);
                nBorderNod = nPeriodic/4; % cause 2D
                Lx  = sum( L(1:nBorderNod) );
                Lxy = sum( L(nBorderNod+1:2*nBorderNod));
                Ly  = sum( L(2*nBorderNod+1 : 3*nBorderNod));
                Ld = L(3*nBorderNod+1 : end); % dirich (2 per + 6 dir)
                switch iBase
                    case 1
                        Lx = Lx + Ld(1) + Ld(2) + Ld(3) + Ld(5);
                        Ly = Ly + Ld(4) + Ld(7);
                    case 2
                        Ly = Ly + Ld(1) + Ld(2) + Ld(4) + Ld(6);
                        Lx = Lx + Ld(3) + Ld(7);
                    case 3
                        Lxy = Lxy + Ld(1) + Ld(2);
                end
                Chomog = [-Lx; -Ly; -Lxy];
           else
               error('Not implemented');
           end
        end

        function n = createFunctionNames(obj, name)
            nStre = numel(obj.uFluc);
            nums = 1:nStre;
            n = cellstr([repmat(name, [nStre,1]), num2str(nums')])';
        end
        
    end

    

end
