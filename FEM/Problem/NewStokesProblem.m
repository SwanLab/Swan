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
        inputBC
        velocityField
        pressureField
    end

    methods (Access = public)

        function obj = NewStokesProblem(cParams)
            obj.init(cParams);
            obj.createVelocityField();
            obj.createPressureField();
            obj.createGeometry();
            obj.createInterpolation();
            obj.createDimensions();
            obj.createDOF();
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

        function createVelocityField(obj) % 1 in old notation
%             bcVelocity = obj.inputBC;
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
%             bcPressure = obj.inputBC;
            bcPressure.dirichlet = obj.inputBC.pressure;
            bcPressure.pointload  = []; % Useless
            s.mesh               = obj.mesh;
            s.ndimf              = 1;
            s.inputBC            = bcPressure;
            s.interpolationOrder = 'LINEAR';
            s.scale              = 'MACRO';
            obj.pressureField = Field(s);
        end

        function createElement(obj)
            obj.element  = Element_Stokes(obj.geometry,obj.mesh,...
                obj.material,obj.dof,obj.problemData,obj.interp,...
                obj.velocityField, obj.pressureField);
        end

        function createSolver(obj)
            velFree  = length(obj.velocityField.boundaryConditions.free);
            prsFree  = length(obj.pressureField.boundaryConditions.free);
            free_dof = [velFree, prsFree];
            s.tol      = 1e-6;
            s.type     = 'Nonlinear';
            s.element  = obj.element;
            s.free_dof = free_dof;
            obj.solver = Solver.create(s);
        end

    end

end
