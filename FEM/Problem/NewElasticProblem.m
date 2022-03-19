classdef NewElasticProblem < handle %NewFEM

%     properties (GetAccess = public, SetAccess = private)
    properties (Access = public)
        variables
    end

    properties (Access = private)
%         material
        nFields
%         interp
        bcApplier
%         quadrature
%         dim
        boundaryConditions
        displacement
    end

    properties (Access = private)
%         mesh
        problemData
        stiffnessMatrix
        stiffnessMatrixRed
        forces
        solver
        geometry
        newBC

        dofsInElem
    end

    properties (Access = protected)
        quadrature
        dim
        material

        vstrain

        mesh, interp % For Homogenization
    end

    methods (Access = public)

        function obj = NewElasticProblem(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.computeDimensions();
            obj.computeDofConnectivity();
            obj.createMaterial();
            obj.computeMaterialProperties();
            obj.createInterpolation();
            obj.createGeometry();
            obj.createBoundaryConditions();
            obj.createNewBoundaryConditions();
            obj.createBCApplier();
            obj.createSolver();
        end

        function solve(obj)
            obj.computeStiffnessMatrix();
            obj.computeForces();
            obj.computeDisplacements();
            obj.computeStrain();
            obj.computeStress();
            obj.computePrincipalDirection();
        end

        function plot(obj)
            s.dim            = obj.dim;
            s.mesh           = obj.mesh;
            s.displacement = obj.variables.d_u;
            plotter = FEMPlotter(s);
            plotter.plot();
        end

        function dim = getDimensions(obj)
            dim = obj.dim;
        end

        function setC(obj, C)
            obj.material.C = C;
        end

        function dvolu = getDvolume(obj)
            dvolu  = obj.mesh.computeDvolume(obj.quadrature);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.nFields = 1;
            obj.mesh        = cParams.mesh;
            pd.scale        = cParams.scale;
            pd.pdim         = cParams.dim;
            pd.ptype        = cParams.type;
            pd.bc.dirichlet = cParams.dirichlet;
            pd.bc.pointload = cParams.pointload;
            if isfield(cParams,'masterSlave')
%                 pd.bc.masterSlave = cParams.masterSlave;
                obj.mesh.computeMasterSlaveNodes();
                pd.bc.masterSlave = obj.mesh.masterSlaveNodes;
            end
            obj.problemData = pd;
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            obj.quadrature = quad;
        end

        function computeDimensions(obj)
            s.ngaus = obj.quadrature.ngaus;
            s.mesh  = obj.mesh;
            s.pdim  = obj.problemData.pdim;
            d       = DimensionVariables(s);
            d.compute();
            obj.dim = d;
        end

        function computeDofConnectivity(obj)
            connec = obj.mesh.connec;
            ndimf  = obj.dim.ndimField;
            nnode  = obj.dim.nnode;
            dofsElem  = zeros(nnode*ndimf,size(connec,1));
            for inode = 1:nnode
                for iunkn = 1:ndimf
                    idofElem   = obj.nod2dof(inode,iunkn);
                    globalNode = connec(:,inode);
                    idofGlobal = obj.nod2dof(globalNode,iunkn);
                    dofsElem(idofElem,:) = idofGlobal;
                end
            end
            obj.dofsInElem = dofsElem;
        end

        function createMaterial(obj)
            s.ptype = obj.problemData.ptype;
            s.pdim  = obj.problemData.pdim;
            s.nelem = obj.mesh.nelem;
            %s.geometry = obj.geometry; % Hyperelastic
            s.mesh  = obj.mesh;
            obj.material = Material.create(s);
        end

        function computeMaterialProperties(obj)
            I = ones(obj.dim.nelem,obj.dim.ngaus);
            s.kappa = .9107*I;
            s.mu    = .3446*I;
            obj.material.compute(s);
        end

        function createInterpolation(obj)
            int = obj.mesh.interpolation;
            obj.interp{1} = int;
        end
       
        function createGeometry(obj)
            q   = obj.quadrature;
            int = obj.interp{1};
            int.computeShapeDeriv(q.posgp);
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(q,int);
            obj.geometry = g;
        end

        function createBoundaryConditions(obj)
            s.dim = obj.dim;
            s.bc  = obj.problemData.bc;
            bc = BoundaryConditions(s);
            bc.compute();
            obj.boundaryConditions = bc;
        end

        function createNewBoundaryConditions(obj)
            s.dim        = obj.dim;
            s.type       = 'Dirichlet';
            s.bc         = obj.problemData.bc;
            s.scale      = obj.problemData.scale;
            s.dofsInElem = obj.dofsInElem;
            bc = NewBoundaryConditions(s);
            bc.compute();
            obj.boundaryConditions = bc;
        end

        function createBCApplier(obj)
            s.BC      = obj.boundaryConditions;
            s.dim     = obj.dim;
            s.nfields = obj.nFields;
            s.scale   = obj.problemData.scale;
%             s.type    = 'Dirichlet'; % defined in Element
            s.type    = ''; % Actually, scale is used...
            obj.bcApplier = BoundaryConditionsApplier.create(s);
        end

        function createSolver(obj)
            obj.solver = Solver.create();
        end

        function computeStiffnessMatrix(obj)
            s.type = 'ElasticStiffnessMatrixOld';
            s.mesh         = obj.mesh;
            s.npnod        = obj.mesh.npnod;
            s.globalConnec = obj.mesh.connec;
            s.dim          = obj.dim;
            s.material     = obj.material;
            LHS = LHSintegrator.create(s);
            K   = LHS.compute();
%             Kred = obj.bcApplier.fullToReducedMatrix(K);
            Kred = obj.boundaryConditions.fullToReducedMatrix(K);
            obj.stiffnessMatrix    = K;
            obj.stiffnessMatrixRed = Kred;
        end

        function computeForces(obj)
            f    = obj.computeExternalForces();
            fRed = obj.reduceForcesMatrix(f);
            obj.forces = fRed;
        end

        function F = computeExternalForces(obj)
            s.dim         = obj.dim;
            s.BC          = obj.boundaryConditions; %dofsInElem, Neumann
            s.dofsInElem  = obj.dofsInElem;
            s.mesh        = obj.mesh;
            s.material    = obj.material;
            s.geometry    = obj.geometry;
            s.dvolume     = obj.getDvolume();
            if isprop(obj, 'vstrain')
                s.vstrain = obj.vstrain;
            end
            fcomp = ForcesComputer(s);
            F = fcomp.compute();
            R = fcomp.computeReactions(obj.stiffnessMatrix);
            obj.variables.fext = F + R;
        end

        function Fred = reduceForcesMatrix(obj, forces)
%             Fred = obj.bcApplier.fullToReducedVector(forces);
            Fred = obj.boundaryConditions.fullToReducedVector(forces);
        end

        function u = computeDisplacements(obj)
            Kred = obj.stiffnessMatrixRed;
            Fred = obj.forces;
            u = obj.solver.solve(Kred,Fred);
%             u = obj.bcApplier.reducedToFullVector(u);
            u = obj.boundaryConditions.reducedToFullVector(u);
            obj.variables.d_u = u;
        end

        function computeStrain(obj)
            s.dim                = obj.dim;
            s.mesh               = obj.mesh;
            s.quadrature         = obj.quadrature;
            s.displacement       = obj.variables.d_u;
            s.interpolation      = obj.interp{1};
            s.dofsInElem         = obj.dofsInElem;
            scomp  = StrainComputer(s);
            strain = scomp.compute();
            obj.variables.strain = strain;
        end

        function computeStress(obj)
            s.C      = obj.material.C;
            s.dim    = obj.dim;
            s.strain = obj.variables.strain;
            scomp  = StressComputer(s);
            stress = scomp.compute();
            obj.variables.stress = stress;
        end

        function computePrincipalDirection(obj)
            stress = obj.variables.stress;
            s.type = obj.problemData.pdim;
            s.eigenValueComputer.type = 'PRECOMPUTED';
            pcomp = PrincipalDirectionComputer.create(s);
            pcomp.compute(stress);
            obj.variables.principalDirections = pcomp.direction;
            obj.variables.principalStress     = pcomp.principalStress;
        end

        function idof = nod2dof(obj, inode, iunkn)
            ndimf = obj.dim.ndimField;
            idof(:,1)= ndimf*(inode - 1) + iunkn;
        end

    end

end