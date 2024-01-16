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
    end

    methods (Access = public)

        function obj = ElasticProblemMicro(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createDisplacementFun();
            obj.createBoundaryConditions();
            obj.createBCApplier();
            obj.createSolver();
        end

        function obj = solve(obj)
            obj.computeStiffnessMatrix();
            obj.computeRHS();
            nCases = obj.material.nstre;
            obj.Chomog = zeros(nCases, nCases);
            for i = 1:nCases
                obj.compDisp(i);
                % obj.computeDisplacements(i);
                obj.computeStrain(i);
                obj.computeStress(i);
                obj.computeFluctuations(i);
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
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            obj.quadrature = quad;
        end

        function createDisplacementFun(obj)
            strdim = regexp(obj.pdim,'\d*','Match');
            nDimf  = str2double(strdim);
            obj.displacementFun = P1Function.create(obj.mesh, nDimf);
        end

        function dim = getFunDims(obj)
            d.ndimf  = obj.displacementFun.ndimf;
            d.nnodes = size(obj.displacementFun.fValues, 1);
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = obj.mesh.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

        function createBoundaryConditions(obj)
            dim = obj.getFunDims();
            bc = obj.inputBC;
            bc.ndimf = dim.ndimf;
            bc.ndofs = dim.ndofs;
            s.mesh  = obj.mesh;
            s.scale = obj.scale;
            s.bc    = {bc};
            s.ndofs = dim.ndofs;
            bc = BoundaryConditions(s);
            bc.compute();
            obj.boundaryConditions = bc;
        end

        function createBCApplier(obj)
            s.mesh = obj.mesh;
            s.boundaryConditions = obj.newBC;
            bc = BCApplier(s);
            obj.BCApplier = bc;
        end

        function createSolver(obj)
            s.type =  'DIRECT';
            obj.solver = Solver.create(s);
        end

        function computeRHS(obj)
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

        function u = compDisp(obj, iVoigt)
            s.type = obj.solverType;
            s.stiffness = obj.stiffness;
            s.forces = obj.forces(:, iVoigt);
            s.boundaryConditions = obj.newBC;
            s.BCApplier = obj.BCApplier;
            pb = ProblemSolver(s); % magic goes here
            u = pb.solve();

            z.mesh    = obj.mesh;
            z.fValues = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            uFeFun = P1Function(z);
            obj.uFun{iVoigt} = uFeFun;

            uSplit = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            obj.displacementFun.fValues = uSplit;
        end

        function u = computeDisplacements(obj, iVoigt)
            [Kred, Fred] = obj.applyBoundaryConditions(iVoigt);
            uRed = obj.solver.solve(Kred,Fred);
            u = obj.boundaryConditions.reducedToFullVector(uRed);

            obj.variables.d_u(:,iVoigt) = u;
            z.mesh    = obj.mesh;
            z.fValues = reshape(u,[obj.mesh.ndim,obj.mesh.nnodes])';
            uFeFun = P1Function(z);
            obj.uFun{iVoigt} = uFeFun;
        end

        function [Kred, Fred] = applyBoundaryConditions(obj,iVoigt)
            bc = obj.boundaryConditions;
            Kred = bc.fullToReducedMatrix(obj.stiffness);
            Fred = bc.fullToReducedVector(obj.forces(:,iVoigt));
        end

        function computeStrain(obj, iVoigt)
            strFun = obj.uFun{iVoigt}.computeSymmetricGradient(obj.quadrature);
            strFun.applyVoigtNotation();
            obj.strainFluctFun{iVoigt} = strFun;
        end

        function computeStress(obj, iVoigt)
            strn  = permute(obj.strainFluctFun{iVoigt}.fValues,[1 3 2]);
            strn2(:,1,:,:) = strn;
            strs =squeeze(pagemtimes(obj.material.C,strn2));
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
            nVoigt  = obj.material.nstre;
            basis   = diag(ones(nVoigt,1));
            vstrain = basis(iVoigt,:);
        end

        function vars = computeFluctuations(obj, iVoigt)
            vstrain = obj.computeVstrain(iVoigt);
            vars  = obj.variables;
            Cmat  = obj.material.C;
            nstre = obj.material.nstre;
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

        function n = createFunctionNames(obj, name)
            nStre = numel(obj.uFun);
            nums = 1:nStre;
            n = cellstr([repmat(name, [nStre,1]), num2str(nums')])';
        end
        
    end

    

end
