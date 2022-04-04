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

    properties (Access = private)
        dim
        boundaryConditions
        M,K, Mr
        epsilon
        solver
        mesh

        dofsInElem
    end

    properties (Access = protected)
        problemData
        bcApplier
    end

    methods (Access = public)
        
        function obj = NewDiffReactProblem(cParams)
            obj.init(cParams);
            obj.createInterpolation();
            obj.computeProblemDimensions();
            obj.computeProblemDofConnectivity();
            obj.createNewBoundaryConditions();
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
%                 bc   = obj.bcApplier;
                bc   = obj.boundaryConditions;
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
                LHS = obj.boundaryConditions.fullToReducedMatrix(LHS);
            end
        end

        function M = getM(obj)
            M = obj.M;
        end

        function M = getK(obj)
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
        
        function computeProblemDimensions(obj)
            m = obj.mesh;
            obj.dim = obj.computeDimensions(m);
        end

        function d = computeDimensions(obj, mesh)
            s.ngaus = [];
            s.mesh  = mesh;
            s.pdim  = obj.problemData.pdim;
            d       = DimensionVariables(s);
            d.compute();
        end

        function computeProblemDofConnectivity(obj)
            d = obj.dim;
            c = obj.mesh.connec;
            obj.dofsInElem = obj.computeDofConnectivity(d, c);
        end

        function dofsElem = computeDofConnectivity(obj, dim, connec)
            ndimf  = dim.ndimField;
            nnode  = dim.nnode;
            dofsElem  = zeros(nnode*ndimf,size(connec,1));
            for inode = 1:nnode
                for iunkn = 1:ndimf
                    idofElem   = obj.nod2dof(inode,iunkn);
                    globalNode = connec(:,inode);
                    idofGlobal = obj.nod2dof(globalNode,iunkn);
                    dofsElem(idofElem,:) = idofGlobal;
                end
            end
        end

        function idof = nod2dof(obj, inode, iunkn)
            ndimf = obj.dim.ndimField;
            idof(:,1)= ndimf*(inode - 1) + iunkn;
        end

        function createNewBoundaryConditions(obj)
            s.dim        = obj.dim;
            s.type       = 'Dirichlet';
            s.bc.dirichlet = [];
            s.bc.pointload = [];
            s.scale      = obj.problemData.scale;
            s.mesh       = obj.mesh;
            s.dofsInElem = obj.dofsInElem;
            bc = NewBoundaryConditions(s);
            bc.compute();
            obj.boundaryConditions = bc;
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
%             switch obj.material
%                 case 'isotropy'
%                     s.type = 'StiffnessMatrix';
%                 case 'anisotropy'
%                     s.type = 'AnisotropicStiffnessMatrix';
%                     for i = 1:size(obj.mesh.connec,1)
%                         s.Celas(:,:,i) = [1, 0; 0, 1]; % Rotation matrix
%                     end
%             end
            s.type = 'StiffnessMatrix';
            s.mesh         = obj.mesh;
            s.npnod        = obj.mesh.npnod;
            s.globalConnec = obj.mesh.connec;
            s.dofsInElem   = obj.dofsInElem;
            s.dim          = obj.dim;
            LHS = LHSintegrator.create(s);
            obj.K = LHS.compute();
        end
        
        function computeMassMatrix(obj)
            s.type         = 'MassMatrix';
            s.quadType     = 'QUADRATICMASS';
            s.mesh         = obj.mesh;
            s.npnod        = obj.mesh.npnod;
            s.dofsInElem   = obj.dofsInElem;
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
                    s.type     = 'MassMatrix';
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
                bfMesh  = boxFaceMesh.mesh;
                gConnec = boxFaceMesh.globalConnec;
                nnode   = size(gConnec,2);
                d.applyNnode(nnode);
                d2 = obj.computeDimensions(bfMesh);
                cParams.mesh = bfMesh;
                cParams.type = 'SIMPLE';
                cParams.globalConnec = gConnec;
                cParams.dofsInElem   = obj.computeDofConnectivity(d2,gConnec);
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