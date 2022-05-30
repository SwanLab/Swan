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
        state
        inputBC
        velocityField
        pressureField
        boundaryConditions
    end

    methods (Access = public)

        function obj = NewStokesProblem(cParams)
            obj.init(cParams);
            obj.createVelocityField();
            obj.createPressureField();
            obj.createBoundaryConditions();
            obj.createGeometry();
            obj.createInterpolation();
            obj.createDOF();
            obj.createElement();
            obj.createSolver();
        end
        
        function computeVariables(obj)
            p.state    = obj.state;
            p.dt         = 0.01; % For transient cases
            p.final_time = 1;    % For transient cases
            x = obj.solver.solve(p);
            obj.variables = obj.element.computeVars(x);
        end

    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.state       = cParams.state;
            pd.scale        = cParams.scale;
            pd.pdim         = cParams.dim;
            pd.ptype        = cParams.type;
            pd.nelem        = cParams.nelem;
            pd.bc.pressure  = cParams.bc.pressure;
            pd.bc.velocity  = cParams.bc.velocity;
            obj.mesh        = cParams.mesh;
            obj.material    = cParams.material;
            obj.problemData = pd;
            obj.fileName    = cParams.fileName;
            obj.inputBC.pressure  = cParams.bc.pressure;
            obj.inputBC.velocity  = cParams.bc.velocity;
            obj.inputBC.pointload = [];
            obj.inputBC.velocityBC = cParams.bc.velocityBC;
            obj.inputBC.forcesFormula = cParams.bc.forcesFormula;
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

        function createDOF(obj)
            obj.dof = DOF_Stokes(obj.fileName,obj.mesh,obj.geometry,obj.interp);
        end

        function createVelocityField(obj) % 1 in old notation
            bcVelocity.dirichlet  = obj.inputBC.velocity;   % Useless
            bcVelocity.pointload  = [];                     % Useless
            bcVelocity.velocityBC = obj.inputBC.velocityBC;
            s.mesh               = obj.mesh;
            s.ndimf              = 2;
            s.inputBC            = bcVelocity;
            s.interpolationOrder = 'QUADRATIC';
            s.scale              = 'MACRO';
            obj.velocityField = Field(s);
        end

        function createPressureField(obj) % 2 in old notation
            bcPressure.dirichlet = obj.inputBC.pressure;
            bcPressure.pointload  = []; % Useless
            s.mesh               = obj.mesh;
            s.ndimf              = 1;
            s.inputBC            = bcPressure;
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'QUADRATIC';
            s.scale              = 'MACRO';
            obj.pressureField = Field(s);
        end

        function createBoundaryConditions(obj)
            vel = obj.velocityField;
            prs = obj.pressureField;
%             bcV = vel.boundaryConditions;
%             bcP = prs.boundaryConditions;
            bcV.dirichlet = vel.translateBoundaryConditions();
            bcV.pointload = [];
            bcV.ndimf     = vel.dim.ndimf;
            bcV.ndofs     = vel.dim.ndofs;
            bcP.dirichlet = obj.inputBC.pressure;
            bcP.pointload = [];
            bcP.ndimf     = prs.dim.ndimf;
            bcP.ndofs     = prs.dim.ndofs;
            ndofs = vel.dim.ndofs + prs.dim.ndofs;
            s.dim   = [];
            s.scale = 'MACRO';
            s.bc    = {bcV, bcP};
            s.ndofs = ndofs; % Stokes
            bc = BoundaryConditions(s);
            obj.boundaryConditions = bc;
        end

        function createElement(obj)
            obj.element  = Element_Stokes(obj.geometry,obj.mesh,...
                obj.material,obj.dof,obj.problemData,obj.interp,...
                obj.velocityField, obj.pressureField, ...
                obj.inputBC.forcesFormula, obj.boundaryConditions);
        end

        function createSolver(obj)
            bc = obj.boundaryConditions;
%             vel = obj.velocityField;
%             prs = obj.pressureField;
%             velCnstr = vel.dim.ndofs - size(vel.inputBC.dirichlet,1);
%             prsCnstr = prs.dim.ndofs - size(prs.inputBC.dirichlet,1);
            free_dof = [length(bc.freeFields{1}), length(bc.freeFields{2})];
            s.tol      = 1e-6;
            s.type     = 'Nonlinear';
            s.element  = obj.element;
            s.free_dof = free_dof;
            obj.solver = Solver.create(s);
        end

    end

end
