classdef ElasticProblemMicroAnisotropic < handle
    
    properties (Access = public)
        uFluc, strain, stress
        Chomog
        Cvoigt
    end

    properties (Access = private)
        mesh
        material
        trialFun
        testFun
        boundaryConditions, bcApplier
        solverType, solverMode, solverCase
        lagrangeMultipliers
        problemSolver
        density
        filter
        fun
    end

    methods (Access = public)

        function obj = ElasticProblemMicroAnisotropic(cParams)
            obj.init(cParams);
            obj.createTrialFun();
            obj.createBCApplier();
            obj.createSolver();
        end

        function obj = solve(obj)
            x = obj.density;
            xD = x.obtainDomainFunction();
            xR = obj.filterField(xD);
            obj.material.setDesignVariable(xR);
            C  = obj.material.obtainTensor();
            dC = obj.material.obtainTensorDerivative();
            %C     = obj.material;
            f     = @(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u)));
            LHS   = IntegrateLHS(f,obj.testFun,obj.trialFun,obj.mesh,'Domain',2);
            for iB = 1:obj.computeNbasis()
                [eB,v] = obj.createDeformationBasis(iB);
                f = @(v) -DDP(SymGrad(v),DDP(C,eB));
                RHS = IntegrateRHS(f,obj.testFun,obj.mesh,'Domain',2);    
                uF{iB}      = obj.computeDisplacement(LHS,RHS,iB);
                strainF{iB} = eB+SymGrad(uF{iB});
                %stressF{iB} = DDP(obj.material, strainF{iB});
                stressF{iB} = DDP(C, strainF{iB});
                ChiB        = Integrator.compute(stressF{iB},obj.mesh,2);
                obj.convertChomogToFourthOrder(ChiB,v,iB);
            end
            obj.uFluc  = uF;
            % Convert tensor to voigt!
            obj.convertChomogToVoigt;
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
            obj.density  = cParams.density;
            obj.mesh     = cParams.mesh;
            obj.material = cParams.material;
            obj.solverType = cParams.solverType;
            obj.solverMode = cParams.solverMode;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.solverCase  = cParams.solverCase;
            obj.filter = cParams.filter;
        end

        function createTrialFun(obj)
            obj.trialFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            obj.testFun = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
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

       function xR = filterField(obj,x)
            nDesVar = length(x);
            xR      = cell(nDesVar,1);
            for i = 1:nDesVar
                xR{i} = obj.filter.compute(x{i},2);
            end
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
            s.solverType = obj.solverType;
            s.solverMode = obj.solverMode;
            s.solver     = obj.solverCase;
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
            uFun.setFValues(full(uSplit));
        end

        function convertChomogToVoigt(obj)
            
            %C_voigt = zeros(3,3);
        
            %map = [1 3; 3 2];  
        %
            %for i = 1:2
            %    for j = 1:2
            %        m = map(i,j);
            %        for k = 1:2
            %            for l = 1:2
            %                n = map(k,l);
            %                C_voigt(m,n) = C_voigt(m,n) + C(i,j,k,l);
            %            end
            %        end
            %    end
            %end
            %obj.Cvoigt = C_voigt;
            
            A = obj.Chomog;
            dim=2;
            if dim == 2
                pairs = [1 1; 2 2; 1 2];
            elseif dim == 3
                pairs = [1 1; 2 2; 3 3; 2 3; 1 3; 1 2];
            else
                error('Tensor must be 2D or 3D.');
            end
            dimVoigt = size(pairs,1);
            
            % ----- Convert 4th-order tensor to Voigt matrix -----
            Avoigt = zeros(dimVoigt);
            for m = 1:dimVoigt
                for n = 1:dimVoigt
                    i = pairs(m,1); j = pairs(m,2);
                    k = pairs(n,1); l = pairs(n,2);
                    Avoigt(m,n) = A(i,j,k,l);
                end
        
            end
            obj.Cvoigt=Avoigt;
        end
        
  
        function convertChomogToFourthOrder(obj,ChiB,v,iB)
            Ch = obj.Chomog;
            v1 = v(iB,1);    v2 = v(iB,2);
            if v1==v2
                Ch(:,:,v1,v2) = ChiB;
            else
                Ch(:,:,v1,v2) = ChiB./2;
                Ch(:,:,v2,v1) = ChiB./2;
            end
            obj.Chomog = Ch;
        end

    end

    

end
