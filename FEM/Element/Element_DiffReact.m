classdef Element_DiffReact < Element
    
    properties
        mesh
        K
        M
        Mr
        epsilon
        interpolation_u
    end
    
    properties (Access = private)
        nstre
        addRobinTerm
        boundaryMesh
    end
    
    methods 
        function obj = Element_DiffReact(mesh,geometry,material,dof,scale, addRobinTerm,bcType,interp,boundaryMesh)
            obj.mesh = mesh;
            obj.addRobinTerm = addRobinTerm;
            obj.bcType = bcType;
            obj.initElement(geometry,mesh,material,dof,scale,interp);
            obj.nstre = 2;
            obj.nfields = 1;
            obj.interpolation_u = interp{1};
            obj.boundaryMesh = boundaryMesh;
            obj.quadrature.computeQuadrature('LINEAR');
            obj.geometry.computeGeometry(obj.quadrature,obj.interpolation_u);
            obj.computeStiffnessMatrix();
            obj.computeMassMatrix();
            obj.computeBoundaryMassMatrix();
        end
        
        function obj = setEpsilon(obj,epsilon)
            obj.epsilon = epsilon;
        end
        
        function LHS = computeLHS(obj)
            if obj.addRobinTerm
                LHS = obj.epsilon^2*obj.K + obj.M + (obj.epsilon)*obj.Mr;
            else
                LHS = obj.epsilon^2*obj.K + obj.M;
                LHS = obj.bcApplier.fullToReducedMatrix(LHS);
            end
        end
        

    end
    
    methods (Access = private)
        
        function computeStiffnessMatrix(obj)
            s.type = 'StiffnessMatrix';
            s.mesh         = obj.mesh;
            s.npnod        = obj.mesh.npnod;
            s.globalConnec = obj.mesh.connec;
            s.dim          = obj.computeDim();
            LHS = LHSintegrator.create(s);
            obj.K = LHS.compute();
        end
        
        function computeMassMatrix(obj)
            s.type         = 'MassMatrix';
            s.quadType     = 'QUADRATICMASS';
            s.mesh         = obj.mesh;
            s.npnod        = obj.mesh.npnod;
            s.globalConnec = obj.mesh.connec;
            s.dim          = obj.computeDim();
            LHS = LHSintegrator.create(s);
            obj.M = LHS.compute();
        end
        
        function computeBoundaryMassMatrix(obj)
            if obj.addRobinTerm
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
            bC      = BoundaryMeshCreator.create(s);
            bMeshes = bC.create();
            bMeshes = obj.boundaryMesh;
            nBoxFaces = numel(bMeshes);
            dim = obj.computeDim();
            for iMesh = 1:nBoxFaces
                boxFaceMesh = bMeshes{iMesh};
                cParams.mesh = boxFaceMesh.mesh;
                cParams.type = 'SIMPLE';
                cParams.globalConnec = boxFaceMesh.globalConnec;
                cParams.npnod        = obj.mesh.npnod;
                cParams.geometryType = obj.mesh.type;
                cParams.dim          = dim;
                params.compositeParams{iMesh} = cParams;
            end
        end
        
        function dim = computeDim(obj)
            s.ngaus = obj.quadrature.ngaus;
            s.mesh  = obj.mesh;
            s.pdim  = '1D';
            dim    = DimensionVariables(s);
            dim.compute();
        end
        
    end
    
    methods(Access = protected) % Only the child sees the function

        function FextSuperficial = computeSuperficialFext(obj)
            FextSuperficial = zeros(obj.nnode*obj.dof.nunkn,1,obj.nelem);
        end
        
        function FextVolumetric = computeVolumetricFext(obj)
            FextVolumetric = zeros(obj.nnode*obj.dof.nunkn,1,obj.nelem);
        end
        
    end

end
