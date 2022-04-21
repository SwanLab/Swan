classdef Hyperelastic_Problem < FEM
    %Hyperelastic_Problem Summary of this class goes here
    % Detailed explanation goes here
    
    %% Public GetAccess properties definition =============================
    properties (GetAccess = public, SetAccess = public)
    end
    
    %% Restricted properties definition ===================================
    properties (Access = private)
        material

        nFields
        interp
        quadrature
    end
    
    %% Public methods definition ==========================================
    methods (Access = public)
        function obj = Hyperelastic_Problem(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();

            obj.createDOF();
            obj.createMaterial();
        end
        
        function preProcess(obj)
            obj.element = Element_Hyperelastic(obj.mesh,obj.geometry,obj.material,obj.dof);
            obj.solver = Solver.create;
        end
        
        function computeVariables(obj)
            tol = 1e-6; % !! This should not be defined in here !!
            for ifield = 1:obj.geometry(1).nfields
                free_dof(ifield) = length(obj.dof.free{ifield});
            end
            x = obj.solve_steady_nonlinear_problem(free_dof,tol);
            obj.variables = obj.element.computeVars(x);
        end

        function init(obj, cParams)
            obj.nFields = 1;
            obj.mesh        = cParams.mesh;
            pd.scale        = cParams.scale;
            pd.pdim         = cParams.dim;
            pd.ptype        = cParams.type;
            pd.bc.dirichlet = cParams.dirichlet;
            pd.bc.pointload = cParams.pointload;
            obj.problemData = pd;
        end
        
        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            obj.quadrature = quad;
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

        function createDOF(obj)
            filename = 'test_hyperelastic';
            mesh = obj.mesh;
            pdim = obj.problemData.pdim;
            nFields = obj.nFields;
            interp = obj.interp;
            obj.dof = DOF_Elastic(filename,mesh,pdim,nFields,interp);
        end

        function createMaterial(obj)
            s.pdim = obj.problemData.pdim;
            s.ptype = obj.problemData.ptype;
            s.nelem = obj.mesh.nelem;
            s.geometry = obj.geometry; % Hyperelastic
            s.mesh  = obj.mesh;
            obj.material = Material.create(s);
        end

    end
end

