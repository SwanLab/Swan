classdef ElasticProblemMicro_Fast < handle
    
    properties (Access = public)
        variables
        uFun, strainFun, stressFun
        strainFluctFun, stressFluctFun
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
        LHS, RHS
    end

    methods (Access = public)

        function obj = ElasticProblemMicro_Fast(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createDisplacementFun();
            obj.createBoundaryConditions();
            obj.createSolver();
            obj.computeStiffnessMatrix();
        end

        function obj = solve(obj)
            obj.computeForces();
            obj.computeDisplacements();
            obj.computeStrain();
            obj.computeStress();
        end

        function obj = computeChomog(obj)
            nVoigt = obj.material.nstre;
            nGaus = obj.quadrature.ngaus;
            nElem = size(obj.material.C,3);
            nPnod = obj.mesh.nnodes; %size(obj.displacementFun.fValues,1);
            strdim = regexp(obj.pdim,'\d*','Match');
            nDimf  = str2double(strdim);
            nDofs = nPnod*nDimf;
            
            tStrn  = zeros(nVoigt,nGaus,nVoigt,nElem);
            tStrss = zeros(nVoigt,nGaus,nVoigt,nElem);
            tDisp  = zeros(nVoigt,nDofs);
            Ch = zeros(nVoigt,nVoigt);
            
            obj.solve();

            for iVoigt = 1:nVoigt
                vstrain = obj.computeVstrain(iVoigt);
                vars = obj.computeStressStrainFluctuations(vstrain, iVoigt);
                Ch(:,iVoigt)         = vars.stress_homog;
                tStrn(iVoigt,:,:,:)  = vars.strain;
                tStrss(iVoigt,:,:,:) = vars.stress;
                tDisp(iVoigt,:)      = obj.variables.d_u(:,iVoigt);
            end
            obj.variables.Chomog  = Ch;
            obj.variables.tstrain = tStrn;
            obj.variables.tstress = tStrss;
            obj.variables.tdisp   = tDisp;
        end

        function computeStiffnessMatrix(obj)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.fun      = obj.displacementFun;
            s.material = obj.material;
            lhs = LHSintegrator.create(s);
            obj.LHS = lhs.compute();
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

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh     = cParams.mesh;
            obj.material = cParams.material;
            obj.scale    = cParams.scale;
            obj.pdim     = cParams.dim;
            obj.inputBC  = cParams.bc;
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

        function createSolver(obj)
            s.type =  'DIRECT';
            obj.solver = Solver.create(s);
        end

        function computeForces(obj)
            s.type = 'ElasticMicroNew';
            s.dim      = obj.getFunDims();
            s.BC       = obj.boundaryConditions;
            s.mesh     = obj.mesh;
            s.material = obj.material;
            s.globalConnec = obj.mesh.connec;
            RHSint = RHSintegrator.create(s);
            rhs = RHSint.compute();
            R = RHSint.computeReactions(obj.LHS);
            obj.variables.fext = rhs + R;
            obj.RHS = rhs;
        end

        function u = computeDisplacements(obj)
            nVoigt = size(obj.RHS,2);
            bc = obj.boundaryConditions;
            Kred = bc.fullToReducedMatrix(obj.LHS);
            Fred = zeros(size(Kred,1), nVoigt);
            for i = 1:nVoigt
                Fred(:,i) = bc.fullToReducedVector(obj.RHS(:,i));
            end
            uRed = obj.solver.solve(Kred,Fred);
            u = zeros(numel(obj.displacementFun.fValues),nVoigt);
            for i = 1:nVoigt
                u(:,i) = bc.reducedToFullVector(uRed(:,i));
                obj.variables.d_u(:,i) = u(:,i);
                z.mesh    = obj.mesh;
                z.fValues = reshape(u(:,i),[obj.mesh.ndim,obj.mesh.nnodes])';
                uFeFun = P1Function(z);
                obj.uFun{i} = uFeFun;
            end
        end

        function computeStrain(obj)
            for i = 1:numel(obj.uFun)
                strFun = obj.uFun{i}.computeSymmetricGradient(obj.quadrature);
                strFun.applyVoigtNotation();
                perm = permute(strFun.fValues, [2 1 3]);
%                 obj.variables.strain = perm;
                obj.strainFluctFun{i} = strFun;
%                 obj.strain = strFun;
            end
        end

        function computeStress(obj)
            for i = 1:numel(obj.strainFluctFun)
                strn  = permute(obj.strainFluctFun{i}.fValues,[1 3 2]);
                strn2(:,1,:,:) = strn;
                strs =squeeze(pagemtimes(obj.material.C,strn2));
                strs = permute(strs, [1 3 2]);
    
                z.mesh       = obj.mesh;
                z.fValues    = strs;
                z.quadrature = obj.quadrature;
                strFun = FGaussDiscontinuousFunction(z);
    
                obj.stress = strFun;
                obj.variables.stress = permute(strFun.fValues, [2 1 3]);
                obj.stressFluctFun{i} = strFun;
            end
        end

        %% 

        function vstrain = computeVstrain(obj, iVoigt)
            nVoigt  = obj.material.nstre;
            basis   = diag(ones(nVoigt,1));
            vstrain = basis(iVoigt,:);
        end

        function prob = solveElasticProblem(obj, prob, vstrain)
            prob.setVstrain(vstrain);
            prob.solve();
        end

        function vars = computeStressStrainFluctuations(obj, vstrain, iVoigt)
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

            vars.stress_fluct = stressFluct;
            vars.strain_fluct = strainFluct;

            vars.stress = stress;
            vars.strain = strain;
            vars.stress_homog = stressHomog;
            obj.variables = vars;
        end
    end

end
