classdef NewElasticProblem < handle %NewFEM

    properties (GetAccess=public, SetAccess = private)
        variables        
    end

    properties (Access = private)
        integrator
        material
        nFields
        interp
        bcApplier

        % Poda 1
        quadrature
        dim   
        LHS
    end

    properties (Access = private)
        fileName
        femData
        mesh
        problemData
        dof
        stiffnessMatrix
        forces
        solver
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
            obj.computeStiffnessMatrix()
            obj.computeForces();
            obj.computeDisplacements();
        end
    
    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.fileName = cParams.fileName;
            obj.nFields = 1;
        end

        function computeStiffnessMatrix(obj)
            obj.LHS = obj.integrator.computeFemLHS(); % lhs need to be adapted
            K = obj.reduceStiffnessMatrix();
            obj.stiffnessMatrix = K;
        end

        function computeForces(obj)
            f    = obj.computeExternalForces();
            fRed = obj.reduceForcesMatrix(f);
            obj.forces = fRed; 
%             R = obj.compute_imposed_displacement_force(obj.K);
%             obj.fext = Fext + R;
%             obj.rhs = obj.integrator.integrate(fNodal);
        end



        function readProblemData(obj)
            obj.createFemData();
            obj.createProblemData();
            obj.mesh = obj.femData.mesh;
        end      

        function createFemData(obj)
            fName       = obj.fileName;
            femReader   = FemInputReader_GiD();
            obj.femData = femReader.read(fName);
        end

        function createProblemData(obj)
            s = obj.femData;
            pd.fileName     = obj.fileName;
            pd.scale        = s.scale;
            pd.pdim         = s.pdim;
            pd.ptype        = s.ptype;
            pd.nelem        = s.mesh.nelem;
            pd.bc.dirichlet = s.dirichlet;
            pd.bc.pointload = s.pointload;
            obj.problemData = pd;
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
            obj.LHS.compute(obj.material.C);
            K = obj.LHS.K;
        end

        % Element_Elastic
        function computeDimensions(obj)
            s.ngaus = obj.quadrature.ngaus;
            s.mesh  = obj.mesh;
            s.pdim  = obj.problemData.pdim;
            d       = DimensionVariables(s);
            d.compute();
            obj.dim = d;
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
            %s.geometry = obj.geometry;
            s.mesh  = obj.mesh;
            obj.material = Material.create(s);
        end

        function createIntegrators(obj)
            s.type         = 'SIMPLE';
            s.mesh         = obj.mesh;
            s.npnod        = obj.mesh.npnod;
            s.globalConnec = obj.mesh.connec;
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
             int = obj.mesh.interpolation;
             %int.computeShapeDeriv(obj.quadrature.posgp);
             obj.interp{1} = int;
         end

        function createSolver(obj)
            obj.solver = Solver.create;
        end
       
        function u = computeDisplacements(obj)
            Kred = obj.stiffnessMatrix;
            Fred = obj.forces;
            u = obj.solver.solve(Kred,Fred);
            u = obj.bcApplier.reducedToFullVector(u);
            obj.variables.d_u = u;
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
                b = zeros(obj.dim.ndof(iField),1);
                nUnkn = obj.dim.nunkn(iField);
                nNode = obj.dim.nnode;
                nDof = obj.dim.ndofPerElement;
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
            FextPoint = zeros(obj.dim.ndof,1);
            if ~isempty(obj.dof.neumann)
                FextPoint(obj.dof.neumann) = obj.dof.neumann_values;
            end
        end

    end

end