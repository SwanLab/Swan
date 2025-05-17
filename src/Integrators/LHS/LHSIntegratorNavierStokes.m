classdef LHSIntegratorNavierStokes < handle

    properties (Access = private)
        mesh
        velocityFun
        velocityField
        material
        LHSS
    end

    methods (Access = public)

        function obj = LHSIntegratorNavierStokes(cParams)
            obj.init(cParams);
        end

        function LHS = compute(obj)
            C   = obj.computeConvectiveMatrix();
            LHS = obj.addCMToLHS(C);
        end

    end

    methods (Access = private)
    
        function init(obj, cParams)
            obj.mesh          = cParams.mesh;
            obj.material      = cParams.material;
            obj.velocityFun   = cParams.velocityFun;
            obj.velocityField = cParams.velocityField;
            obj.LHSS          = cParams.LHSS;
        end

        function c = computeConvectiveMatrix(obj)
            s.type  = 'Convective';
            s.mesh  = obj.mesh;
            s.test  = obj.velocityFun;
            s.trial = obj.velocityFun;
            s.material = obj.material;
            C = LHSIntegrator.create(s);
            c = C.compute(obj.velocityField);
        end

        function LHS = addCMToLHS(obj,C)
            LHS  = obj.LHSS;
            nVel = size(C, 1);
            LHS(1:nVel, 1:nVel) = LHS(1:nVel, 1:nVel) + C;
        end

    end

end