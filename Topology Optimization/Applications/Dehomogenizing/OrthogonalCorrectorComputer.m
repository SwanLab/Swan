classdef OrthogonalCorrectorComputer < handle

    properties (Access = private)
        shiftingFunction
        orthogonalCorrector
    end

    properties (Access = private)
        mesh
        interpolator
        corrector
    end

    methods (Access = public)

        function obj = OrthogonalCorrectorComputer(cParams)
            obj.init(cParams)
        end

        function oC = compute(obj)
            obj.createShifting();
            obj.createOrthogonalCorrector();
            oC = obj.orthogonalCorrector;
        end

        function plot(obj)
            obj.corrector.plot();            
            obj.shiftingFunction.plot();
            obj.orthogonalCorrector.plot();
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.interpolator = cParams.interpolator;
            obj.corrector    = cParams.corrector;
        end

        function createShifting(obj)
            s.mesh     = obj.mesh;
            s.fValue   = obj.corrector.fValues;
            s.rhsType = 'ShapeDerivative';
            s.interpolator = obj.interpolator;
            m = MinimumDiscGradFieldWithVectorInH1(s);
            f = m.solve();
            s.mesh = obj.mesh;
            s.fValues(1,:,:) = f';
            sh = P1DiscontinuousFunction(s);
            obj.shiftingFunction = sh;
        end

        function createOrthogonalCorrector(obj)
            phi = obj.corrector.fValues;
            fD  = obj.shiftingFunction.fValues;
            phi = phi - fD;
            s.mesh    = obj.mesh;
            s.fValues = phi;
            cF = P1DiscontinuousFunction(s);  
            obj.orthogonalCorrector = cF;         
        end

    end

end