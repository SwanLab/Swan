classdef NewElasticProblem < NewFEM

    properties (Access = private)
        integrator
        material
        nFields
        interp
        bcApplier

        % Poda 1
        quadrature
        dim
    end

    methods (Access = public)

        function obj = NewElasticProblem(cParams)
            obj.init(cParams);
            obj.readProblemData();
            obj.createQuadrature();
            obj.computeDimensions();
            obj.createMaterial();
            obj.createInterpolation();
            obj.createBCApplier();
            obj.createIntegrators();
            obj.createSolver();
        end
        
        function computeVariables(obj)
            obj.integrator.computeLHS(); % lhs need to be adapted
            Kred = obj.reduceStiffnessMatrix();
            forces = obj.computeExternalForces();
%             R = obj.compute_imposed_displacement_force(obj.K);
%             obj.fext = Fext + R;
            Fred = obj.reduceForcesMatrix(forces);
%             obj.rhs = obj.integrator.integrate(fNodal);
            u = obj.solver.solve(Kred,Fred);
            obj.variables = obj.processVars(u);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.fileName = cParams.fileName;
            obj.nFields = 1;
        end

        function Fred = reduceForcesMatrix(obj, forces)
            Fred = obj.bcApplier.fullToReducedVector(forces);
        end

        function Kred = reduceStiffnessMatrix(obj)
            K = obj.computeStiffnessMatrixSYM();
            Kred = obj.bcApplier.fullToReducedMatrix(K);
        end

        function K = computeStiffnessMatrixSYM(obj)
            obj.computeC();
            obj.integrator.StiffnessMatrix.compute(obj.material.C);
            K = obj.integrator.StiffnessMatrix.K;
        end

        % Element_Elastic
        function dim = computeDimensions(obj)
            dim                = DimensionVariables();
            dim.nnode          = obj.mesh.nnode;
            dim.nunkn          = obj.createNUnkn();
            dim.nstre          = obj.createNstre();
            dim.ndof           = obj.mesh.npnod*dim.nunkn;
            dim.nelem          = obj.mesh.nelem;
            dim.ndofPerElement = dim.nnode*dim.nunkn;
            dim.ngaus          = obj.quadrature.ngaus;
            dim.nentries       = dim.nelem*(dim.ndofPerElement)^2;
            dim.ndim           = obj.createNdim();
            obj.dim = dim;
        end

        % IsotropicElasticMaterial
        function computeC(obj)
            pdim = obj.problemData.pdim;
            switch pdim
                case '2D'
                    obj.computeC2D();
                case '3D'
                    obj.computeC3D();
            end
        end
        
        function computeC2D(obj)
            I = ones(obj.mesh.nelem,obj.quadrature.ngaus);
            kappa = .9107*I;
            mu    = .3446*I;
            nElem = size(mu,1);
            nGaus = size(mu,2);
            m = mu;
            l = kappa - mu;
            C = zeros(obj.dim.nstre,obj.dim.nstre,nElem,nGaus);
            C(1,1,:,:)= 2*m+l;
            C(1,2,:,:)= l;
            C(2,1,:,:)= l;
            C(2,2,:,:)= 2*m+l;
            C(3,3,:,:)= m;
            obj.material.C = C;
        end
        
        function computeC3D(obj)
            I = ones(obj.mesh.nelem,obj.quadrature.ngaus);
            kappa = .9107*I;
            mu    = .3446*I;
            nElem = size(mu,1);
            nGaus = size(mu,2);
            m = mu;
            l = kappa - 2/3*mu;
            C = zeros(obj.dim.nstre,obj.dim.nstre,nElem,nGaus);
            C(1,1,:,:) = 2*m + l;
            C(2,2,:,:) = 2*m + l;
            C(3,3,:,:) = 2*m + l;
            C(1,2,:,:) = l;
            C(2,1,:,:) = l;
            C(1,3,:,:) = l;
            C(3,1,:,:) = l;
            C(3,2,:,:) = l;
            C(2,3,:,:) = l;
            C(4,4,:,:) = m;
            C(5,5,:,:) = m;
            C(6,6,:,:) = m;
            obj.material.C = C;
        end

        function createMaterial(obj)
            s.ptype = obj.problemData.ptype;
            s.pdim  = obj.problemData.pdim;
            s.nelem = obj.mesh.nelem;
            s.geometry = obj.geometry;
            s.mesh  = obj.mesh;
            obj.material = Material.create(s);
        end

        function createIntegrators(obj)
            s.type = 'SIMPLE';
            s.mesh  = obj.mesh;
            s.npnod = obj.mesh.npnod;
            s.fileName     = obj.fileName;
            s.globalConnec = obj.mesh.connec;
            s.problemData  = obj.problemData;
            s.bcApplier    = obj.bcApplier;
            s.dim          = obj.dim;
            obj.integrator = Integrator.create(s);
        end

        function createBCApplier(obj)
            obj.dof = DOF_Elastic(obj.fileName,obj.mesh,obj.problemData.pdim,obj.nFields,obj.interp);
            cParams.nfields = obj.nFields;
            cParams.dof = obj.dof;
            cParams.scale = obj.problemData.scale;
            cParams.type = 'Dirichlet'; % defined in Element
            obj.bcApplier = BoundaryConditionsApplier.create(cParams);
        end
        
        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            obj.quadrature = quad;
        end

         function createInterpolation(obj)
             int = Interpolation.create(obj.mesh,'LINEAR');
             %int.computeShapeDeriv(obj.quadrature.posgp);
             obj.interp{1} = int;
         end

        function createSolver(obj)
            obj.solver = Solver.create;
        end

        function variables = processVars(obj, uL)
            variables.d_u = obj.computeDisplacements(uL);
        end

        function u = computeDisplacements(obj, usol)
            u = obj.bcApplier.reducedToFullVector(usol);
        end

        function Fext = computeExternalForces(obj)
            FextSuperficial = obj.computeSuperficialFext;
            FextVolumetric  = obj.computeVolumetricFext;
            FextSupVol = {FextSuperficial + FextVolumetric};
            FextSupVol = obj.AssembleVector(FextSupVol);
            FextPoint = obj.computePunctualFext();
            Fext = FextSupVol +  FextPoint;
        end

        function FextSuperficial = computeSuperficialFext(obj)
            d = obj.dim;
            nnode = d.nnode;
            nunkn = d.nunkn;
            nelem = obj.mesh.nelem;
            FextSuperficial = zeros(nnode*nunkn,1,nelem);
        end
        
        function FextVolumetric = computeVolumetricFext(obj)
            d = obj.dim;
            nnode = d.nnode;
            nunkn = d.nunkn;
            nelem = obj.mesh.nelem;
            FextVolumetric = zeros(nnode*nunkn,1,nelem);
        end

        function b = AssembleVector(obj,b_elem_cell)
            nfields = 1;
            for iField = 1:nfields
                bElem = b_elem_cell{iField,1};
                b = zeros(obj.dof.ndof(iField),1);
                nUnkn = obj.dof.nunkn(iField);
                nNode = obj.interp{iField}.nnode;
                nDof = nNode*nUnkn;
                nGaus = size(bElem,2);
                for iDof = 1:nDof
                    for igaus = 1:nGaus
                        c = squeeze(bElem(iDof,igaus,:));
                        idof_elem = obj.dof.in_elem{iField}(iDof,:);
                        b = b + sparse(idof_elem,1,c',obj.dof.ndof(iField),1);
                    end
                end
                b_global{iField,1} = b;
            end
            b=cell2mat(b_global);
        end

        function FextPoint = computePunctualFext(obj)
            %Compute Global Puntual Forces (Not well-posed in FEM)
            FextPoint = zeros(obj.dof.ndof,1);
            if ~isempty(obj.dof.neumann)
                FextPoint(obj.dof.neumann) = obj.dof.neumann_values;
            end
        end

        function nUnkn = createNUnkn(obj)
            pdim = obj.problemData.pdim;
            switch pdim
                case '2D'
                    nUnkn = 2;
                case '3D'
                    nUnkn = 3;
            end
        end

        function ndim = createNdim(obj)
            pdim = obj.problemData.pdim;
            switch pdim
                case '2D'
                    ndim = 2;
                case '3D'
                    ndim = 3;
            end
        end

        function nstre = createNstre(obj)
            pdim = obj.problemData.pdim;
            switch pdim
                case '2D'
                    nstre = 3;
                case '3D'
                    nstre = 6;
            end
        end

    end

end