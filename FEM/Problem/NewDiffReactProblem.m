classdef NewDiffReactProblem < handle %FEM
    
    properties (Access = protected)
        isRobinTermAdded
        bcApplierType
        interp
        boundaryMesh
    end
    
    properties (GetAccess = public, SetAccess = protected) %FEM
        geometry
        variables
    end

    % new properties
    properties (Access = private)
        dim
        boundaryConditions
        M,K, Mr
        epsilon
%         bcApplier
        solver
%         problemData
        mesh
    end

    properties (Access = protected)
        problemData
        bcApplier
    end

    methods (Access = public)
        
        function obj = NewDiffReactProblem(cParams)
            obj.init(cParams);
            obj.createInterpolation();
            obj.computeDimensions();
            obj.createBoundaryConditions();
            obj.createBCApplier();
            obj.createGeometry();
            obj.createSolver();
            obj.computeStiffnessMatrix();
            obj.computeMassMatrix();
            obj.computeBoundaryMassMatrix();
        end
        
        function computeVariables(obj,x)
            if obj.isRobinTermAdded
                LHS  = obj.computeLHS();
                xReg = obj.solver.solve(LHS,x);
                obj.variables.x = xReg;
            else
                bc   = obj.bcApplier;
                xRed = bc.fullToReducedVector(x);
                LHS  = obj.computeLHS();
                xReg = obj.solver.solve(LHS,xRed);
                obj.variables.x = bc.reducedToFullVector(xReg);
            end
        end
        

        function obj = setEpsilon(obj,epsilon)
            obj.epsilon = epsilon;
        end
        
        function LHS = computeLHS(obj)
            if obj.isRobinTermAdded
                LHS = obj.epsilon^2*obj.K + obj.M + (obj.epsilon)*obj.Mr;
            else
                LHS = obj.epsilon^2*obj.K + obj.M;
                LHS = obj.bcApplier.fullToReducedMatrix(LHS);
            end
        end

        function M = getM(obj) %new
            M = obj.M;
        end

        function M = getK(obj) %new
            M = obj.K;
        end

    end
    
    methods (Access = protected)
        
        function setScale(obj)
            obj.problemData.scale = 'MACRO';
        end
        
    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.setupRobinTermAndBCApplier(cParams)
            obj.setupFromMesh(cParams);
            obj.problemData.pdim = '1D';
            obj.setScale();
        end

        function createInterpolation(obj)
            int = obj.mesh.interpolation;
            obj.interp{1} = int;
        end
        
        function computeDimensions(obj)
            s.ngaus = [];
            s.mesh  = obj.mesh;
            s.pdim  = obj.problemData.pdim;
            d       = DimensionVariables(s);
            d.compute();
            obj.dim = d;
        end

        function createBoundaryConditions(obj)
            s.dim          = obj.dim;
            s.globalConnec = obj.mesh.connec;
            s.bc.dirichlet = [];
            s.bc.pointload = [];
            bc = BoundaryConditions(s);
            bc.compute();
            obj.boundaryConditions = bc;
        end

        function createBCApplier(obj)
            s.BC      = obj.boundaryConditions;
            s.dim     = obj.dim;
            s.scale   = obj.problemData.scale;
            s.type    = obj.bcApplierType;
            s.nfields = 1;
            s.mesh = obj.mesh;
            obj.bcApplier = BoundaryConditionsApplier.create(s);
        end

        function createGeometry(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(q,obj.interp{1});
            obj.geometry = g;
        end
        
        function createSolver(obj)
            obj.solver = Solver.create();
        end

        function computeStiffnessMatrix(obj)
            s.type = 'StiffnessMatrix';
            s.mesh         = obj.mesh;
            s.npnod        = obj.mesh.npnod;
            s.globalConnec = obj.mesh.connec;
            s.dim          = obj.dim;
            LHS = LHSintegrator.create(s);
            obj.K = LHS.compute();
        end
        
        function computeMassMatrix(obj)
            s.type         = 'MassMatrix';
            s.quadType     = 'QUADRATICMASS';
            s.mesh         = obj.mesh;
            s.npnod        = obj.mesh.npnod;
            s.globalConnec = obj.mesh.connec;
            s.dim          = obj.dim;
            LHS = LHSintegrator.create(s);
            obj.M = LHS.compute();
        end
        
        function computeBoundaryMassMatrix(obj)
            if obj.isRobinTermAdded
                cParams = obj.createIntegratorParams();
                nInt = numel(cParams.compositeParams);
                ndof = cParams.compositeParams{1}.dim.ndof;
                LHS = sparse(ndof,ndof);
                for iInt = 1:nInt
                    s = cParams.compositeParams{iInt};
                    s.type = 'MassMatrix';
                    s.quadType = 'LINEAR';
                    lhs = LHSintegrator.create(s);
                    LHSadd = lhs.compute();
                    LHS = LHS + LHSadd;
                end
                obj.Mr = LHS;
            end
        end
        
        function params = createIntegratorParams(obj)
            params.type  = 'COMPOSITE';
            params.npnod = obj.mesh.npnod;
            s.backgroundMesh = obj.mesh;
            s.dimension      = 1:s.backgroundMesh.ndim;
            s.type           = 'FromReactangularBox';
            bC = BoundaryMeshCreator.create(s);
            bC.create();
            bMeshes = obj.boundaryMesh;
            nBoxFaces = numel(bMeshes);
            d = obj.dim;
            for iMesh = 1:nBoxFaces
                boxFaceMesh = bMeshes{iMesh};
                cParams.mesh = boxFaceMesh.mesh;
                cParams.type = 'SIMPLE';
                cParams.globalConnec = boxFaceMesh.globalConnec;
                cParams.npnod        = obj.mesh.npnod;
                cParams.geometryType = obj.mesh.type;
                cParams.dim          = d;
                params.compositeParams{iMesh} = cParams;
            end
        end

        function setupRobinTermAndBCApplier(obj, cParams)
            if isfield(cParams,'isRobinTermAdded')
                obj.isRobinTermAdded = cParams.isRobinTermAdded;
                obj.bcApplierType = cParams.bcApplierType;
            else
                obj.isRobinTermAdded = false;
                obj.bcApplierType = '';
            end
        end

        function setupFromMesh(obj,s)
            obj.mesh = s.mesh;
            if isfield(s,'fileName') % UnffitedIntegration
                obj.problemData.fileName = s.fileName;
                obj.createBoundaryMesh(s.fileName);
            end
        end
        
        function createBoundaryMesh(obj,fileName)
            run(fileName);
            if exist('External_border_nodes','var') && ~isempty(External_border_nodes)
                s.borderNodes    = External_border_nodes;
                s.borderElements = External_border_elements;
                s.backgroundMesh = obj.mesh;
                s.type = 'FromData';
                b = BoundaryMeshCreator.create(s);
                obj.boundaryMesh = b.create();
            else
                s.backgroundMesh = obj.mesh;
                s.dimension = 1:obj.mesh.ndim;
                s.type = 'FromReactangularBox';
                bC = BoundaryMeshCreator.create(s);
                obj.boundaryMesh = bC.create();
            end
        end
    
    end
    
end