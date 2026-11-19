classdef ImportedSpecimenGeometry < SpecimenGeometry

    methods (Access = public)

        function obj = ImportedSpecimenGeometry(cParams)
            obj.init(cParams);
            obj.setGripsFromBoundingBox();
            obj.setSideGrips(cParams.fixedSideHeight,cParams.nCells);
        end

    end

    methods (Static, Access = public)

        function cP = curvedBeamParams(nCells)
            here = fileparts(mfilename('fullpath'));
            cP.fileName   = fullfile(here,'curvedBeam_boundary.csv');
            cP.rotation   = 90;
            cP.hingeWidth = 1.2;
            cP.fixedSideHeight = 15;
            cP.nCells     = nCells;
        end

    end

    methods (Access = protected)

        function points = computeRefinePoints(obj)
            points = obj.shape.Vertices(all(isfinite(obj.shape.Vertices),2),:);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            p = readmatrix(cParams.fileName);
            unitCell = rotate(polyshape(p(:,1),p(:,2)),cParams.rotation,[0 0]);
            obj.shape           = obj.tile(unitCell,cParams.nCells);
            obj.hingeWidth      = cParams.hingeWidth;
            obj.cornerLigament  = cParams.hingeWidth;
            obj.hingeTips       = zeros(0,2);
            obj.latticeVertices = zeros(0,2);
            [~,yLim] = boundingbox(unitCell);
            obj.sideLength      = diff(yLim);
        end

        function tiled = tile(~,unitCell,nCells)
            [xLim,yLim] = boundingbox(unitCell);
            cells = cell(nCells^2,1);
            for i = 0:nCells-1
                for j = 0:nCells-1
                    v = unitCell.Vertices + [i*diff(xLim), j*diff(yLim)];
                    cells{i*nCells+j+1} = polyshape(round(v,9));
                end
            end
            tiled = union([cells{:}]);
        end

    end

end
