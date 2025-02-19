classdef ElasticProblemMicro < handle
    
    properties (Access = public)
        variables
        uFun, strainFun, stressFun
        strainFluctFun, stressFluctFun
        Chomog
    end

    properties (Access = private)
        mesh
        pdim
        material
        quadrature
        displacementFun
        boundaryConditions, bcApplier
        strain, stress
        stiffness, forces

        solverType, solverMode, solverCase
        lagrangeMultipliers

        problemSolver
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
            oX     = zeros(obj.getDimensions().ndimf,1);
            nCases = size(obj.material.evaluate(oX),1);
            obj.Chomog = zeros(nCases, nCases);
            for i = 1:nCases
                obj.computeDisplacement(i);
                obj.computeStrain(i);
                obj.computeStress(i);
                obj.computeChomogContribution(i);
            end
        end

        function computeStiffnessMatrix(obj)
            ndimf = obj.displacementFun.ndimf;
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.test     = LagrangianFunction.create(obj.mesh,ndimf, 'P1');
            s.trial    = obj.displacementFun;
            s.material = obj.material;
            s.quadratureOrder = 2;
            lhs = LHSIntegrator.create(s);
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
%             interp.computeShapeDeriv(obj.quadrature.posgp);
        end

        function quad = getQuadrature(obj)
            quad = obj.quadrature;
        end

        function updateMaterial(obj, mat)
            obj.material = mat;
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
            obj.pdim     = cParams.dim;
            obj.solverType = cParams.solverType;
            obj.solverMode = cParams.solverMode;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.solverCase  = cParams.solverCase;
        end

        function createQuadrature(obj)
            quad = Quadrature.create(obj.mesh, 1);
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

        function computeForces(obj)
            s.fun  = obj.displacementFun;
            s.type = 'ElasticMicro';
            s.dim      = obj.getFunDims();
            s.BC       = obj.boundaryConditions;
            s.mesh     = obj.mesh;
            s.material = obj.material;
            s.globalConnec = obj.mesh.connec;
            RHSint = RHSIntegrator.create(s);
            rhs = RHSint.compute();
            R = RHSint.computeReactions(obj.stiffness);
            obj.variables.fext = rhs + R;
            obj.forces = rhs;
        end

        function u = computeDisplacement(obj, iVoigt)
            s.stiffness = obj.stiffness;
            s.forces    = obj.forces(:, iVoigt);
            s.iVoigt    = iVoigt;
            s.nVoigt    = size(obj.forces,2);
            [u, L]      = obj.problemSolver.solve(s);
            obj.lagrangeMultipliers = L;
            z.mesh    = obj.mesh;
            z.fValues = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            z.order   = 'P1';
            uFeFun = LagrangianFunction(z);
            obj.uFun{iVoigt} = uFeFun;

            uSplit = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            obj.displacementFun.setFValues(uSplit);
        end

        function computeStrain(obj, iVoigt)
            nCases    = size(obj.Chomog,1);
            e         = zeros(nCases,1,1);
            e(iVoigt) = 1;
            strn      = SymGrad(obj.uFun{iVoigt});

            obj.strainFluctFun{iVoigt} = strn;
            s.operation                = @(xV) e+strn.evaluate(xV);
            s.mesh                     = obj.mesh;
            s.ndimf                    = strn.ndimf;
            obj.strainFun{iVoigt}      = DomainFunction(s);
        end

        function computeStress(obj, iVoigt)
            obj.stress = DDP(obj.material, obj.strainFluctFun{iVoigt});
%             obj.variables.stress = permute(strFun.fValues, [2 1 3]);
            obj.stressFluctFun{iVoigt} = DDP(obj.material, obj.strainFluctFun{iVoigt});
        end

        %% 

        function vstrain = computeVstrain(obj, iVoigt)
            oX      = zeros(obj.getDimensions().ndimf,1);
            nVoigt  = size(obj.material.evaluate(oX),1);
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
                oX    = zeros(obj.getDimensions().ndimf,1);
                nstre = size(obj.material.evaluate(oX),1);
                nelem = size(Cmat,4);
                ngaus = obj.quadrature.ngaus;
                dV = obj.mesh.computeDvolume(obj.quadrature)';
                strainFluct = obj.strainFluctFun{iVoigt}.evaluate(xV);
                stressFluct = obj.stressFluctFun{iVoigt}.evaluate(xV);
                
                stress = zeros(nstre,ngaus,nelem);
                strain = zeros(nstre,ngaus,nelem);
                stressHomog = zeros(nstre,1);
                
                for igaus = 1:ngaus
                    strain(1:nstre,igaus,:) = vstrain'.*ones(nstre,1,nelem) + strainFluct(1:nstre,ngaus,:);
                    for istre = 1:nstre
                        for jstre = 1:nstre
                            Cij  = squeeze(Cmat(istre,jstre,igaus,:));
                            C    = squeeze(Cij);
                            strs = squeeze(stress(istre,igaus,:));
                            strn = squeeze(strain(jstre,igaus,:));
                            stress(istre,igaus,:) = strs + C.* strn;
                        end
                        strs = squeeze(stress(istre,igaus,:));
                        stressHomog(istre) = stressHomog(istre) + (strs)'*dV(:,igaus);
                    end
                end
    
                obj.Chomog(:,iVoigt) = stressHomog;
                
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
