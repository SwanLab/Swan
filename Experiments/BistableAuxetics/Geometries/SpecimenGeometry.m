classdef SpecimenGeometry < handle

    properties (GetAccess = public, SetAccess = protected)
        shape
        hingeTips
        bottomGrips
        topGrips
        sideGrips
        height
        sideLength
        hingeWidth
        latticeVertices
        cornerLigament
        center
    end

    methods (Access = public)

        function ls = computeLevelSetFunction(obj,mesh)
            ls = LagrangianFunction.create(mesh,1,'P1');
            ls.setFValues(PolyshapeSignedDistance(obj.shape).compute(mesh.coord));
        end

        function mesh = createMesh(obj,refinement,workDir,name)
            s = obj.computeMeshSizes(refinement);
            s.boundaryPoints = obj.latticeVertices;
            s.refinePoints   = obj.computeRefinePoints();
            s.isBackground   = false;
            mesh = obj.runMesher(s,workDir,name);
        end

        function mesh = createBackgroundMesh(obj,refinement,workDir,name)
            s = obj.computeBackgroundMeshSizes(refinement);
            s.boundaryPoints = zeros(0,2);
            s.refinePoints   = zeros(0,2);
            s.isBackground   = true;
            mesh = obj.runMesher(s,workDir,name);
        end

        function rho = computeDensity(obj,mesh)
            chi = obj.createCharacteristicFunction(mesh);
            s.filterType = 'LUMP';
            s.mesh       = mesh;
            s.trial      = LagrangianFunction.create(mesh,1,'P1');
            rho = Filter.create(s).compute(chi,3);
        end

        function s = computeBarrierParameters(obj)
            s.activationDistance = 0.2*obj.hingeWidth;
        end

    end

    methods (Abstract, Access = protected)
        points = computeRefinePoints(obj)
    end

    methods (Access = private)

        function mesh = runMesher(obj,s,workDir,name)
            s.shape           = obj.shape;
            s.refinementShape = obj.computeRefinementShape();
            s.workDir         = workDir;
            s.name            = name;
            mesh = GmshPolyshapeMesher(s).compute();
        end

        function chi = createCharacteristicFunction(obj,mesh)
            s.fun      = obj.computeLevelSetFunction(mesh);
            s.mesh     = mesh;
            s.type     = 'LevelSet';
            s.plotting = false;
            chi = DesignVariable.create(s).obtainDomainFunction();
            chi = chi{1};
        end

    end

    methods (Access = protected)

        function s = computeMeshSizes(obj,refinement)
            w = obj.hingeWidth;
            L = obj.sideLength;
            s.sizeMin = refinement*w/4;
            s.sizeMax = refinement*L/10;
            s.distMin = 2*w;
            s.distMax = L/3;
        end

        function w = computeCutWidth(obj)
            w = obj.hingeWidth;
        end

        function shape = computeRefinementShape(obj)
            shape = obj.shape;
        end

        function s = computeBackgroundMeshSizes(obj,refinement)
            h = refinement*min(obj.computeCutWidth(),obj.hingeWidth)/4;
            s.sizeMin = h;
            s.sizeMax = 8*h;
            s.distMin = obj.computeCutWidth();
            s.distMax = obj.computeCutWidth() + 8*h;
        end

        function setGripsFromBoundingBox(obj)
            [xLim,yLim] = boundingbox(obj.shape);
            margin = obj.sideLength;
            obj.bottomGrips = [xLim(1)-margin, xLim(2)+margin, yLim(1)];
            obj.topGrips    = [xLim(1)-margin, xLim(2)+margin, yLim(2)];
            obj.height      = diff(yLim);
            obj.center      = [mean(xLim),mean(yLim)];
        end

        function setSideGrips(obj,sideHeight,nRows)
            [xLim,yLim] = boundingbox(obj.shape);
            yStart = yLim(1) + obj.sideLength*(0:nRows-1)';
            sides  = [yStart, yStart+sideHeight];
            obj.sideGrips = [sides, xLim(1)*ones(nRows,1); sides, xLim(2)*ones(nRows,1)];
        end

    end

end
