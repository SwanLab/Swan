classdef OrthogonalCorrectorComputer < handle

    properties (Access = public)

    end

    properties (Access = private)

        shiftingValue
        orthogonalCorrectorValue
    end

    properties (Access = private)
        mesh
        interpolator
        correctorValue
    end

    methods (Access = public)

        function obj = OrthogonalCorrectorComputer(cParams)
            obj.init(cParams)
        end

        function cF = compute(obj)
            obj.createShifting();
            obj.createOrthogonalCorrector();
            c(1,:,:) = obj.orthogonalCorrectorValue';
            s.mesh    = obj.mesh;
            s.fValues = c;
            cF = P1DiscontinuousFunction(s);
        end

        function plot(obj)
            obj.plotFieldDG((obj.orthogonalCorrectorValue))
            obj.plotFieldDG((obj.shiftingValue))

            figure()
            m = obj.mesh.createDiscontinousMesh();
            x = m.coord(:,1);
            y = m.coord(:,2);
            z = abs(obj.shiftingValue');
            %figure()
            tricontour(m.connec,x,y,z,linspace(min(z(:)),max(z(:)),30))
            %   view(0,90)
            %    colorbar
            %    shading interp
        end

        function plotFieldDG(obj,f)
            figure()
            s.mesh  = obj.mesh.createDiscontinousMesh();
            s.field = transpose(f);
            n = NodalFieldPlotter(s);
            n.plot();
            shading interp
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.interpolator       = cParams.interpolator;
            obj.correctorValue     = cParams.correctorValue;
        end

        function createShifting(obj)
            s.mesh     = obj.mesh;
            s.fValue   = obj.correctorValue;
            s.rhsType = 'ShapeDerivative';
            s.interpolator = obj.interpolator;
            m = MinimumDiscGradFieldWithVectorInH1(s);
            f = m.solve();
            obj.shiftingValue = f;
        end

        function createOrthogonalCorrector(obj)
            phi = obj.correctorValue;
            fD  = obj.shiftingValue;
            phi = phi - fD;
            obj.orthogonalCorrectorValue = phi;
        end


    end

end