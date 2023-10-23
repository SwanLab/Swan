classdef LHSintegrator_DiffReactRobin < handle

    properties (GetAccess = public, SetAccess = private)
        test, trial
        M
        K
        Mr
    end

    properties (Access = private)
        mesh
    end

    methods (Access = public)

        function obj = LHSintegrator_DiffReactRobin(cParams)
            obj.init(cParams);
            obj.computeStiffnessMatrix(cParams);
            obj.computeMassMatrix();
            obj.computeBoundaryMassMatrix();
        end

        function LHS = compute(obj, epsilon)
            LHS = epsilon^2*obj.K + obj.M + epsilon*obj.Mr;
        end

    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
            obj.test  = cParams.trial;
            obj.trial = cParams.trial;
        end

        function computeStiffnessMatrix(obj,cParams)
            s      = cParams; % For anisotropic stiffness
            s.test  = obj.test;
            s.trial  = obj.trial;
            s.type = cParams.stiffType;
            s.mesh = obj.mesh;
            LHS    = LHSintegrator.create(s);
            obj.K  = LHS.compute();
        end

        function computeMassMatrix(obj)
            s.type = 'MassMatrix';
            s.mesh = obj.mesh;
            s.test  = obj.test;
            s.trial = obj.trial;
            s.quadratureOrder = 'QUADRATICMASS';
            LHS     = LHSintegrator.create(s);
            obj.M   = LHS.compute();
        end

        function computeBoundaryMassMatrix(obj)
%             g.mesh    = obj.mesh;
%             g.fValues = zeros(size(obj.mesh.coord,1),1);
%             f = P1Function(g);
            s.type  = 'BoundaryMassMatrix';
            s.mesh  = obj.mesh;
%             s.fun   = f;
            s.quadratureOrder = 'QUADRATICMASS';
            LHS     = LHSintegrator.create(s);
            obj.Mr  = LHS.compute();
        end

    end

end