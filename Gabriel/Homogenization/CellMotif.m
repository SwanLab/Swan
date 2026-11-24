classdef CellMotif < handle

    properties (Access = private)
        type
        pnorm
        center
        tolerance
        maxIterations
        scaleBracket
        aspectRatio
    end

    methods (Access = public)

        function obj = CellMotif(cParams)
            obj.type          = cParams.type;
            obj.pnorm         = cParams.pnorm;
            obj.center        = cParams.center(:);
            obj.tolerance     = cParams.tolerance;
            obj.maxIterations = cParams.maxIterations;
            obj.scaleBracket  = cParams.scaleBracket;
            obj.aspectRatio   = 1;
            if isfield(cParams,'aspectRatio')
                obj.aspectRatio = cParams.aspectRatio;
            end
        end

        function dens = computeDensity(obj,mesh,T,scale)
            if scale == 0
                dens = obj.createConstantDensity(mesh,1);
            elseif isinf(scale)
                dens = obj.createConstantDensity(mesh,0);
            else
                dens = obj.createLevelSetDensity(mesh,T,scale);
            end
        end

        function phi = evaluateLevelSet(obj,x,T,scale)
            nPoints = size(x,2);
            if scale == 0
                phi = ones(1,nPoints);
            elseif isinf(scale)
                phi = -ones(1,nPoints);
            else
                g   = GeometricalFunction(obj.createGeometryParameters(scale));
                fH  = obj.createPeriodicHandle(g.getHandle(),T);
                phi = reshape(fH(reshape(x,size(x,1),nPoints,1)),1,nPoints);
            end
        end

        function rho = computeVolumeFraction(~,mesh,dens)
            rho = Integrator.compute(dens,mesh,2)/mesh.computeVolume();
        end

        function scale = computeScale(obj,referenceMesh,rho)
            if rho >= 1
                scale = 0;
            elseif rho <= 0
                scale = Inf;
            else
                scale = obj.bisectScale(referenceMesh,rho);
            end
        end

    end

    methods (Access = private)

        function scale = bisectScale(obj,mesh,rho)
            T        = eye(mesh.ndim);
            scaleMin = obj.scaleBracket(1);
            scaleMax = obj.scaleBracket(2);
            while obj.computeRho(mesh,T,scaleMax) > rho
                scaleMax = 2*scaleMax;
            end
            for i = 1:obj.maxIterations
                scale      = 0.5*(scaleMin+scaleMax);
                rhoCurrent = obj.computeRho(mesh,T,scale);
                if abs(rhoCurrent-rho) < obj.tolerance
                    return
                end
                if rhoCurrent > rho
                    scaleMin = scale;
                else
                    scaleMax = scale;
                end
            end
            scale = 0.5*(scaleMin+scaleMax);
        end

        function rho = computeRho(obj,mesh,T,scale)
            dens = obj.computeDensity(mesh,T,scale);
            rho  = obj.computeVolumeFraction(mesh,dens);
        end

        function dens = createConstantDensity(~,mesh,value)
            dens = LagrangianFunction.create(mesh,1,'P1');
            dens.setFValues(value*ones(mesh.nnodes,1));
        end

        function dens = createLevelSetDensity(obj,mesh,T,scale)
            ls = obj.computeLevelSet(mesh,T,scale);
            sU.backgroundMesh = mesh;
            sU.boundaryMesh   = mesh.createBoundaryMesh;
            uMesh = UnfittedMesh(sU);
            uMesh.compute(ls);
            chi = CharacteristicFunction.create(uMesh);
            sF.trial = LagrangianFunction.create(mesh,1,'P1');
            sF.mesh  = mesh;
            filter = FilterLump(sF);
            dens = filter.compute(chi,2);
        end

        function ls = computeLevelSet(obj,mesh,T,scale)
            g         = GeometricalFunction(obj.createGeometryParameters(scale));
            s.fHandle = obj.createPeriodicHandle(g.getHandle(),T);
            s.ndimf   = 1;
            s.mesh    = mesh;
            aFun      = AnalyticalFunction(s);
            phi       = aFun.project('P1');
            ls        = -phi.fValues;
        end

        function fH = createPeriodicHandle(obj,fRef,T)
            Tinv = T\eye(size(T));
            c    = obj.center;
            wrap = @(xi) c + mod(xi - c + 0.5,1) - 0.5;
            fH   = @(x) fRef(wrap(pagemtimes(Tinv,x)));
        end

        function gPar = createGeometryParameters(obj,scale)
            gPar.type        = obj.type;
            gPar.pnorm       = obj.pnorm;
            gPar.xCoorCenter = obj.center(1);
            gPar.yCoorCenter = obj.center(2);
            switch obj.type
                case 'Circle'
                    gPar.radius = scale/2;
                case 'Square'
                    gPar.length = scale;
                case 'Rectangle'
                    gPar.xSide = scale*obj.aspectRatio;
                    gPar.ySide = scale/obj.aspectRatio;
                case 'Ellipse'
                    gPar.type  = 'SmoothRectangle';
                    gPar.xSide = scale;
                    gPar.ySide = scale;
                    gPar.pnorm = 2;
                case 'SmoothHexagon'
                    gPar.radius = scale;
                    gPar.normal = [0 1;sqrt(3)/2 1/2;sqrt(3)/2 -1/2];
                otherwise
                    error('CellMotif:UnknownType','Unknown motif type: %s',obj.type);
            end
        end

    end

end