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
            nTerms = 3*nBasis;
            % obj.Chomog = zeros(nBasis, nBasis);

            for iB = 1:nTerms

                if iB <= 3
                    s = obj.createDeformationBasis(iB);
                    iBasis = iB;
                elseif iB <= 6
                    iBasis = iB - 3;
                    coord = 1;
                    s = obj.createSndOrderDeformationBasis(coord, iBasis);
                else
                    iBasis = iB - 6;
                    coord = 2;
                    s = obj.createSndOrderDeformationBasis(coord, iBasis);
                end
                
                strainB     = s;
                RHS         = obj.computeRHS(strainB,LHS);
                uF{iB}      = obj.computeDisplacement(LHS,RHS,iBasis,nBasis);
                strainF{iB} = strainB+SymGrad(uF{iB});
                stressF{iB} = DDP(obj.material, strainF{iB});
                %  Ch(:,iB)    = obj.computeChomog(stressF{iB},iB);
                uM     = obj.computeTotal(strainB);
                uT{iB} = uF{iB} + uM;
            end

            obj.uFluc  = uF;
            obj.strain = strainF;
            obj.stress = stressF;
            % obj.Chomog = Ch;
        end

        function uT = computeTotal(obj,strainB)
            Y  = AnalyticalFunction.create(@(x) x,2,obj.mesh);
            uMF = @(xV) obj.obtainMacroscopicDisplacement(xV,strainB,Y);
            uT = AnalyticalFunction.create(uMF,2,obj.mesh);
        end 

        function uM = obtainMacroscopicDisplacement(xV,strain,Y)
            e  = strain.evaluate(xV);            
            y = Y.evaluate(xV);
            y1 = y(1,:,:);
            y2 = y(2,:,:);
            ex = e(1,:,:);
            ey = e(2,:,:);
            exy = e(3,:,:);
            uM(1,:,:) = ex.*y1 + exy.*y2;
            uM(2,:,:) = exy.*y1 + ey.*y2;
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

        function s = createDeformationBasis(obj,iBasis)
            nBasis = obj.computeNbasis();
            sV = zeros(nBasis,1);
            sV(iBasis) = 1;
            s = ConstantFunction.create(sV,obj.mesh);
        end

        function s = createSndOrderDeformationBasis(obj,coord,iBasis) 
            nBasis = obj.computeNbasis();
            sV = zeros(nBasis,1);
            sV(iBasis) = 1;
            eta = ConstantFunction.create(sV,obj.mesh);
            Y1  = AnalyticalFunction.create(@(x) x(coord,:,:),1,obj.mesh);
            s = eta.*Y1;
            % s.ndimf = nBasis;
        end

        function nBasis = computeNbasis(obj)
            nDim = obj.mesh.ndim;
            nBasis=nDim*(nDim+1)/2;
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

        function Chomog = computeChomog(obj,stress,iBase)
            if strcmp(obj.solverMode, 'DISP')
                Chomog = computeChomogFromLagrangeMultipliers(obj,iBase);
            else
                Chomog = Integrator.compute(stress,obj.mesh,2);
            end
        end

        function Chomog = computeChomogFromLagrangeMultipliers(obj,iBase)
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
        end

    end

    

end
