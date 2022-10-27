classdef FunElasticProblem < handle
    
    % Just a concept. Does not work.
    properties (Access = public)
        
    end

    properties (Access = private) % Calculated
        LHS, RHS, solver
        displacement
    end

    properties (Access = private) % Inputs
        pdim
        mesh
        scale
        material
        boundaryConditions
    end

    methods (Access = public)

        function obj = FunElasticProblem(cParams)
            obj.init(cParams);
            obj.createDisplacementFunction();
%             obj.createBoundaryConditions();
            obj.createSolver();
        end

        function solve(obj)
            obj.computeStiffnessMatrix();
%             obj.computeForces();
%             obj.computeDisplacements();
%             obj.computeStrain();
%             obj.computeStress();
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.pdim     = cParams.dim;
            obj.mesh     = cParams.mesh;
            obj.scale    = cParams.scale;
%             obj.inputBC  = cParams.bc;
            obj.material = cParams.material;
        end

        function createDisplacementFunction(obj)
            strdim = regexp(obj.pdim,'\d*','Match');
            nDimf  = str2double(strdim);
            nNodes = size(obj.mesh.coord,1);
            s.ndimf   = nDimf;
            s.connec  = obj.mesh.connec;
            s.type    = obj.mesh.type;
            s.fValues = zeros(nNodes,nDimf);
            f = P1Function(s);
            obj.displacement = f;
        end

        function createBoundaryConditions(obj)
            s.boundaryConditions = obj.boundaryConditions;
            s.mesh  = obj.mesh;
            s.scale = obj.scale;
            bc = BoundaryConditions(s);
            bc.compute();
            obj.boundaryConditions = bc;
        end

        function createSolver(obj)
            s.type =  'DIRECT';
            obj.solver = Solver.create(s);
        end

        function computeStiffnessMatrix(obj)
            s.type     = 'FunElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.fun      = obj.displacement;
            s.material = obj.material;
            lhs = LHSintegrator.create(s);
            obj.LHS = lhs.compute();
        end

        function computeForces(obj)
            s.type = 'Elastic';
            s.scale    = obj.scale;
            s.dim      = obj.displacementField.dim;
            s.BC       = obj.boundaryConditions;
            s.mesh     = obj.mesh;
            s.material = obj.material;
            s.globalConnec = obj.displacementField.connec;
            if isprop(obj, 'vstrain')
                s.vstrain = obj.vstrain;
            end
            RHSint = RHSintegrator.create(s);
            rhs = RHSint.compute();
            R = RHSint.computeReactions(obj.LHS);
            obj.variables.fext = rhs + R;
            obj.RHS = rhs;
        end

        function u = computeDisplacements(obj)
            bc = obj.boundaryConditions;
            Kred = bc.fullToReducedMatrix(obj.LHS);
            Fred = bc.fullToReducedVector(obj.RHS);
            u = obj.solver.solve(Kred,Fred);
            u = bc.reducedToFullVector(u);
            obj.displacement.fValues = u;
        end

        function computeStrain(obj)
            obj.strain = obj.displacement.computeGradient();
        end

        function computeStress(obj)
            obj.stress = pagemtimes(obj.material.C, obj.strain);
        end
        

    end

end
