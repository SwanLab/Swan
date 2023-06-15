classdef ThermalProblem < handle
    
    properties (Access = public)
        variables
    end

    properties (Access = private)
        temperature
        boundaryConditions
        problemData
        stiffnessMatrix
        stiffnessMatrixRed
        forces
        solver
        geometry
        quadrature        
    end

    properties (Access = private)
        mesh        
        scale
        pdim
        inputBC          
    end

    methods (Access = public)

        function obj = ThermalProblem(cParams)
            obj.init(cParams);
            obj.createTemperatureFun();
            obj.createBoundaryConditions();
            obj.createSolver();
        end

        function solve(obj)
            obj.computeStiffnessMatrix();
            obj.computeForces();
            obj.computeTemperatures();
        end

        function dvolu = getDvolume(obj)
            dvolu  = obj.mesh.computeDvolume(obj.quadrature);
        end

        function print(obj,fileName)
            s = obj.createPostProcessDataBase(fileName);
            s.pdim = '2D';
            s.name = 'temp';
            s.ptype     = obj.problemData.ptype;
            s.ndim      = 2;obj.dim.ndim;
            %s.pdim      = obj.problemData.pdim;
            postprocess = Postprocess('ScalarNodal',s);
            q = obj.quadrature;
            d.fields = obj.variables.d_u;
            d.quad = q;
            postprocess.print(1,d);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh        = cParams.mesh;
            obj.scale       = cParams.scale;
            obj.pdim        = cParams.dim;
            obj.inputBC     = cParams.bc;
        end

        function createTemperatureFun(obj)
            strdim = regexp(obj.pdim,'\d*','Match');
            nDimf  = str2double(strdim);
            t = P1Function.create(obj.mesh, nDimf);
            obj.temperature = t;
        end

        function createBoundaryConditions(obj)
            dim = obj.getFunDims();
            bc = obj.inputBC;
            bc.ndimf = dim.ndimf;
            bc.ndofs = dim.ndofs;
            s.mesh  = obj.mesh;
            s.scale = obj.scale;
            s.bc    = {bc};
            s.ndofs = dim.ndofs;
            bc = BoundaryConditions(s);
            bc.compute();
            obj.boundaryConditions = bc;
        end  

        function dim = getFunDims(obj)
            d.ndimf  = obj.temperature.ndimf;
            d.nnodes = size(obj.temperature.fValues, 1);
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = obj.mesh.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

        function createSolver(obj)
            s.type = 'DIRECT';
            obj.solver = Solver.create(s);
        end

        

        function computeStiffnessMatrix(obj)
            s.type     = 'StiffnessMatrix';
            s.mesh     = obj.mesh;
            s.fun      = obj.temperature;
            s.material = obj.material;
            lhs = LHSintegrator.create(s);
            obj.LHS = lhs.compute();
        end

        function computeForces(obj)
            f    = obj.computeExternalForces();
            fRed = obj.reduceForcesMatrix(f);
            obj.forces = fRed;
        end

        function F = computeExternalForces(obj)
            s.dim         = obj.dim;
            s.BC          = obj.boundaryConditions; % Neumann
            s.mesh        = obj.mesh;
            s.material    = [];
            s.geometry    = obj.geometry;
            s.dvolume     = obj.getDvolume();
            s.globalConnec = obj.mesh.connec;
            fcomp = ForcesComputer(s);
            F = fcomp.compute();
            R = fcomp.computeReactions(obj.stiffnessMatrix);
            obj.variables.fext = F + R;
        end

        function Fred = reduceForcesMatrix(obj, forces)
            Fred = obj.boundaryConditions.fullToReducedVector(forces);
        end

        function t = computeTemperatures(obj)
            Kred = obj.stiffnessMatrixRed;
            Fred = obj.forces;
            t = obj.solver.solve(Kred,Fred);
            t = obj.boundaryConditions.reducedToFullVector(t);
            obj.variables.d_u = t;
        end

        function d = createPostProcessDataBase(obj,fileName)
            dI.mesh    = obj.mesh;
            dI.outFileName = fileName;
            dI.pdim = '2D';
            dI.ptype = 'THERMAL';
            ps = PostProcessDataBaseCreator(dI);
            d = ps.create();
        end
        
    end

end