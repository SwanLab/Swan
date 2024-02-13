classdef ElasticProblemMicro < handle
    
    properties (Access = public)
        variables
        uFun, strainFun, stressFun
        strainFluctFun, stressFluctFun
        Chomog
    end

    properties (Access = private)
        mesh
        scale
        pdim
        inputBC
        material
        quadrature
        displacementFun
        solver
        boundaryConditions
        strain, stress
        stiffness, forces

        solverType, solverMode
        newBC, BCApplier
        lagrangeMultipliers
    end

    methods (Access = public)

        function obj = ElasticProblemMicro(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createDisplacementFun();
            obj.createBCApplier();
            obj.createSolver();
        end

        function obj = solve(obj)
            obj.computeStiffnessMatrix();
            obj.computeForces();
            nCases = size(obj.material.evaluate([0;0]),1);
            obj.Chomog = zeros(nCases, nCases);
            for i = 1:nCases
                obj.computeDisplacement(i);
                obj.computeStrain(i);
                obj.computeStress(i);
                obj.computeChomogContribution(i);
            end
        end

        function computeStiffnessMatrix(obj)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.fun      = obj.displacementFun;
            s.material = obj.material;
            lhs = LHSintegrator.create(s);
            obj.stiffness = lhs.compute();
        end

        function v = computeGeometricalVolume(obj)
            v = 1;%sum(sum(obj.geometry.dvolu));
        end

        function dvolu = getDvolume(obj)
            dvolu  = obj.mesh.computeDvolume(obj.quadrature);
        end

        function setMatProps(obj,s)
           obj.material.compute(s);
        end

        function mesh = getMesh(obj)
            mesh  = obj.mesh;
        end
        
        function interp = getInterpolation(obj)
            interp  = obj.mesh.interpolation;
            interp.computeShapeDeriv(obj.quadrature.posgp);
        end

        function quad = getQuadrature(obj)
            quad = obj.quadrature;
        end

        function setC(obj, C)
            obj.material.C = C;
        end

        function dim = getDimensions(obj)
            strdim = regexp(obj.pdim,'\d*','Match');
            nDimf  = str2double(strdim);
            d.ndimf  = nDimf;
            d.nnodes = obj.mesh.nnodes;
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = obj.mesh.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

        function [fun, funNames] = getFunsToPlot(obj)
            if isempty(obj.stressFun)
                fun = [];
                funNames = [];
            else
                dispN = obj.createFunctionNames('displacement');
                strsN = obj.createFunctionNames('stress');
                strnN = obj.createFunctionNames('strain');
                fun = {obj.uFun{:}, obj.strainFun{:}, obj.stressFun{:}};
                funNames = {dispN{:}, strsN{:}, strnN{:}};
                obj.uFun = {};
                obj.strainFun = {};
                obj.stressFun = {};
            end
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh     = cParams.mesh;
            obj.material = cParams.material;
            obj.scale    = cParams.scale;
            obj.pdim     = cParams.dim;
            obj.inputBC  = cParams.bc;
            obj.solverType  = cParams.solverType;
            obj.solverMode  = cParams.solverMode;
            obj.newBC = cParams.newBC;
            obj.boundaryConditions = cParams.boundaryConditions;
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            obj.quadrature = quad;
        end

        function createDisplacementFun(obj)
            strdim = regexp(obj.pdim,'\d*','Match');
            nDimf  = str2double(strdim);
            obj.displacementFun = LagrangianFunction.create(obj.mesh, nDimf, 'P1');
        end

        function dim = getFunDims(obj)
            d.ndimf  = obj.displacementFun.ndimf;
            d.nnodes = size(obj.displacementFun.fValues, 1);
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = obj.mesh.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

        function createBCApplier(obj)
            s.mesh = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            bc = BCApplier(s);
            obj.BCApplier = bc;
        end

        function createSolver(obj)
            s.type =  'DIRECT';
            obj.solver = Solver.create(s);
        end

        function computeForces(obj)
            s.fun  = obj.displacementFun;
            s.type = 'ElasticMicro';
            s.dim      = obj.getFunDims();
            s.BC       = obj.boundaryConditions;
            s.mesh     = obj.mesh;
            s.material = obj.material;
            s.globalConnec = obj.mesh.connec;
            RHSint = RHSintegrator.create(s);
            rhs = RHSint.compute();
            R = RHSint.computeReactions(obj.stiffness);
            obj.variables.fext = rhs + R;
            obj.forces = rhs;
        end

        function u = computeDisplacement(obj, iVoigt)
            s.solverType = obj.solverType;
            s.solverMode = obj.solverMode;
            s.stiffness = obj.stiffness;
            s.forces = obj.forces(:, iVoigt);
            s.boundaryConditions = obj.boundaryConditions;
            s.boundaryConditions.iVoigt = iVoigt;
            s.boundaryConditions.nVoigt = size(obj.forces,2);
            s.BCApplier = obj.BCApplier;
            pb = ProblemSolver(s); % magic goes here
            [u, L] = pb.solve();

            obj.lagrangeMultipliers = L;
            z.mesh    = obj.mesh;
            z.fValues = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            z.order   = 'P1';
            uFeFun = LagrangianFunction(z);
            obj.uFun{iVoigt} = uFeFun;

            uSplit = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            obj.displacementFun.fValues = uSplit;
        end

        function computeStrain(obj, iVoigt)
            strFun = obj.uFun{iVoigt}.computeSymmetricGradient(obj.quadrature);
            obj.strainFluctFun{iVoigt} = strFun.obtainVoigtFormat();
        end

        function computeStress(obj, iVoigt)
            xV = obj.quadrature.posgp;
            Cmat = obj.material.evaluate(xV);
            strn  = permute(obj.strainFluctFun{iVoigt}.fValues,[1 3 2]);
            strn2(:,1,:,:) = strn;
            strs =squeeze(pagemtimes(Cmat,strn2));
            strs = permute(strs, [1 3 2]);

            z.mesh       = obj.mesh;
            z.fValues    = strs;
            z.quadrature = obj.quadrature;
            strFun = FGaussDiscontinuousFunction(z);

            obj.stress = strFun;
            obj.variables.stress = permute(strFun.fValues, [2 1 3]);
            obj.stressFluctFun{iVoigt} = strFun;
        end

        %% 

        function vstrain = computeVstrain(obj, iVoigt)
            nVoigt  = size(obj.material.evaluate([0;0]),1);
            basis   = diag(ones(nVoigt,1));
            vstrain = basis(iVoigt,:);
        end

        function vars = computeChomogContribution(obj, iVoigt)
            if strcmp(obj.solverMode, 'DISP')
                L = obj.lagrangeMultipliers;
                nPeriodic = length(obj.boundaryConditions.periodic_leader);
                nBorderNod = nPeriodic/4; % cause 2D
                Lx  = sum( L(1:nBorderNod) );
                Lxy = sum( L(nBorderNod+1:2*nBorderNod));
                Ly  = sum( L(2*nBorderNod+1 : 3*nBorderNod));
                Ld = L(3*nBorderNod+1 : end); % dirich (2 per + 6 dir)
                switch iVoigt
                    case 1
                        Lx = Lx + Ld(1) + Ld(2) + Ld(3) + Ld(5);
                        Ly = Ly + Ld(4) + Ld(7);
                    case 2
                        Ly = Ly + Ld(1) + Ld(2) + Ld(4) + Ld(6);
                        Lx = Lx + Ld(3) + Ld(7);
                    case 3
                        Lxy = Lxy + Ld(1) + Ld(2);
                end
                obj.Chomog(iVoigt,:) = [-Lx; -Ly; -Lxy];

            else
                vstrain = obj.computeVstrain(iVoigt);
                vars  = obj.variables;
                xV    = obj.quadrature.posgp;
                Cmat  = obj.material.evaluate(xV);
                nstre = size(obj.material.evaluate([0;0]),1);
                nelem = size(Cmat,3);
                ngaus = obj.quadrature.ngaus;
                dV = obj.mesh.computeDvolume(obj.quadrature)';
                strainFluct = permute(obj.strainFluctFun{iVoigt}.fValues, [2 1 3]);
                stressFluct = permute(obj.stressFluctFun{iVoigt}.fValues, [2 1 3]);
                
                stress = zeros(ngaus,nstre,nelem);
                strain = zeros(ngaus,nstre,nelem);
                stressHomog = zeros(nstre,1);
                
                for igaus = 1:ngaus
                    strain(igaus,1:nstre,:) = vstrain.*ones(1,nstre,nelem) + strainFluct(igaus,1:nstre,:);
                    for istre = 1:nstre
                        for jstre = 1:nstre
                            Cij  = squeeze(Cmat(istre,jstre,:,igaus));
                            C    = squeeze(Cij);
                            strs = squeeze(stress(igaus,istre,:));
                            strn = squeeze(strain(igaus,jstre,:));
                            stress(igaus,istre,:) = strs + C.* strn;
                        end
                        strs = squeeze(stress(igaus,istre,:));
                        stressHomog(istre) = stressHomog(istre) + (strs)'*dV(:,igaus);
                    end
                end
    
                obj.Chomog(:,iVoigt) = stressHomog;
    
                a.mesh       = obj.mesh;
                a.fValues    = permute(stress, [2 1 3]);
                a.quadrature = obj.quadrature;
                obj.stressFun{iVoigt} = FGaussDiscontinuousFunction(a);
    
                a.mesh       = obj.mesh;
                a.fValues    = permute(strain, [2 1 3]);
                a.quadrature = obj.quadrature;
                obj.strainFun{iVoigt} = FGaussDiscontinuousFunction(a);
    
    
                vars.stress_fluct = stressFluct;
                vars.strain_fluct = strainFluct;
    
                vars.stress = stress;
                vars.strain = strain;
                vars.stress_homog = stressHomog;
                obj.variables = vars;
            end
        end

        function n = createFunctionNames(obj, name)
            nStre = numel(obj.uFun);
            nums = 1:nStre;
            n = cellstr([repmat(name, [nStre,1]), num2str(nums')])';
        end
        
    end

    

end
