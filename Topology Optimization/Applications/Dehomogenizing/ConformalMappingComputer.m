classdef ConformalMappingComputer < handle


    properties (Access = private)
        phi
        interpolator
        phiMapping
        totalCorrector
    end

    properties (Access = private)
        orientation
        orientationVector
        mesh
        dilation
    end

    methods (Access = public)

        function obj = ConformalMappingComputer(cParams)
            obj.init(cParams);
            obj.computeOrientationVector();
            obj.createInterpolator();
        end

        function phiV = compute(obj)
            obj.computeMappingWithSingularities();
            phiV = obj.phi;
        end

        function plot(obj)
            for iDim = 1:obj.mesh.ndim
                obj.phi{iDim}.plotContour();
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.orientation = cParams.theta;
            obj.mesh        = cParams.mesh;
            obj.dilation    = cParams.dilation;
        end

        function createInterpolator(obj)
            s.mesh        = obj.mesh;
            s.orientation = [cos(obj.orientation),sin(obj.orientation)];
            s = SymmetricContMapCondition(s);
            sC = s.computeCondition();
            obj.interpolator = sC;
        end

        function computeMappings(obj)
            nDim = obj.mesh.ndim;
            for iDim = 1:nDim
                bI    = obj.orientationVector{iDim}.fValues;
                phiD  = obj.computeMapping(bI);
                phiV(iDim,:,:) = phiD.fValues;
            end
            s.fValues = phiV;
            s.mesh    = obj.mesh;
            phiF = P1DiscontinuousFunction(s);
            obj.phiMapping = phiF;
        end

        function computeTotalCorrector(obj)
           s.mesh = obj.mesh;
           s.orientationVector = obj.orientationVector;
           s.interpolator = obj.interpolator;
           s.phiMapping   = obj.phiMapping;
           tC = TotalCorrectorComputer(s);           
           obj.totalCorrector = tC.compute();
        end


        function computeMappingWithSingularities(obj)
            obj.computeMappings();
            obj.computeTotalCorrector();
            psiTs = P1DiscontinuousFunction.sum(obj.phiMapping,obj.totalCorrector);
            obj.phi = psiTs;
        end

        function phi = computeMapping(obj,fValues)
            s.fValues = fValues';
            s.mesh    = obj.mesh;
            s.rhsType = 'ShapeDerivative';
            s.interpolator = obj.interpolator;
            problem   = MinimumDiscGradFieldWithVectorInL2(s);
            phi = problem.solve();
        end


        function computeOrientationVector(obj)
            er = exp(obj.dilation);
            erCos = er.*cos(obj.orientation);
            erSin = er.*sin(obj.orientation);
            b1(:,1) = erCos;
            b1(:,2) = erSin;
            b2(:,1) = -erSin;
            b2(:,2) = erCos;
            b(:,:,1) = b1;
            b(:,:,2) = b2;
            for iDim = 1:obj.mesh.ndim
                s.fValues = b(:,:,iDim);
                s.mesh   = obj.mesh;
                bf = P1Function(s);
                obj.orientationVector{iDim} = bf;
            end
        end

    end

end
