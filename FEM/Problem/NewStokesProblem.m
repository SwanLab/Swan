classdef NewStokesProblem < handle

    properties (Access = public)
        variables
    end
    
    properties (Access = private)
        interp
        problemData
        geometry
        mesh
        dof
        element
        material
        solver
        fileName

        dim
        boundaryConditions
    end

    methods (Access = public)

        function obj = NewStokesProblem(cParams)
            obj.init(cParams);
            obj.createGeometry();
            obj.createInterpolation();
            obj.createDimensions();
            obj.createDOF();
            obj.createBoundaryConditions();
            obj.createMaterial();
            obj.createElement();
            obj.createSolver();
        end
        
        function computeVariables(obj)
            p.state    = 'Steady';
%             p.state      = 'Transient';
%             p.dt         = 0.01; % For transient cases
%             p.final_time = 1;    % For transient cases
            x = obj.solver.solve(p);
            obj.variables = obj.element.computeVars(x);
        end

    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.fileName = cParams.fileName;
%             obj.problemData.fileName = cParams.fileName;
            obj.problemData.scale = cParams.scale;
            obj.problemData.pdim  = cParams.dim;
            obj.problemData.ptype = cParams.type;
            obj.problemData.nelem = cParams.nelem;
            obj.problemData.bc.dirichlet = cParams.bc.dirichlet;
            obj.problemData.bc.pointload = cParams.bc.pointload;
            obj.mesh = cParams.mesh;
        end

        function createGeometry(obj)
            s.mesh = obj.mesh;
            obj.geometry    = Geometry.create(s);
            obj.geometry(2) = Geometry.create(s);
        end
        
        function createInterpolation(obj)
            interpU = 'QUADRATIC';
            interpP = 'LINEAR';
            obj.interp{1}=Interpolation.create(obj.mesh,interpU);
            obj.interp{2}=Interpolation.create(obj.mesh,interpP);
        end

        function createDimensions(obj)
            v.type = 'Vector';
            v.fieldName = 'v';
            v.mesh = obj.mesh;
            v.ndimf = 2;
            vDim = DimensionVariables.create(v);
            vDim.compute();
            p.type = 'Scalar';
            p.name = 'p';
            p.mesh = obj.mesh;
            p.ndimf = 1;
            pDim = DimensionVariables.create(p);
            obj.dim = {vDim, pDim};
            % This is wrong. The mesh creates a linear interpolation by
            % default, with no option to change it.
        end

        function createDOF(obj)
            obj.dof = DOF_Stokes(obj.fileName,obj.mesh,obj.geometry,obj.interp);
        end

        function createBoundaryConditions(obj)
            PP = Preprocess;
            pD = obj.problemData;
            [fixnodes,forces,~,~] = PP.getBC_fluids...
                (obj.fileName, obj.mesh, obj.geometry,obj.interp);
            % fixnodes = dirichlet_data;
            % forces = neumann_data;
            obj.createBCvelocity();
            obj.createBCpressure();
        end

        function createBCvelocity(obj, dirichlet, neumann)
            s.dim        = obj.dim{1};
            s.mesh       = obj.mesh; % nope
            s.scale      = obj.problemData.scale;
            s.bc         = obj.problemData.bc;
            bcV = BoundaryConditions(s);
            bcV.compute();
            obj.boundaryConditions{1} = bcV;
        end

        function createBCpressure(obj)
            s.dim = obj.dim{2};
            s.mesh   = obj.mesh;
            s.scale  = obj.problemData.scale;
            s.bc     = obj.problemData.bc;
            bcP = BoundaryConditions(s);
            bcP.compute();
            obj.boundaryConditions{2} = bcP;
        end

        function createMaterial(obj)
            cParams.nelem = obj.mesh.nelem;
            mat = Material_Stokes(cParams);
            mat.compute();
            obj.material = mat;

        end

        function createElement(obj)
            obj.element  = Element_Stokes(obj.geometry,obj.mesh,obj.material,obj.dof,obj.problemData,obj.interp, obj.dim);
        end

        function createSolver(obj)
            nFields = numel(obj.interp);
            for ifield = 1:nFields
                free_dof(ifield) = length(obj.dof.free{ifield});
            end
            s.tol      = 1e-6;
            s.type     = 'Nonlinear';
            s.element  = obj.element;
            s.free_dof = free_dof;
            obj.solver = Solver.create(s);
        end

    end

end
