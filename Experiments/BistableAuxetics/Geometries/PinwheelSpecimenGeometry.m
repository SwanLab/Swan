classdef PinwheelSpecimenGeometry < SpecimenGeometry

    properties (Access = private)
        nSides
        cuts
        cutRoots
        t
        theta
        s
        nCells
        tiling
    end

    methods (Access = public)

        function obj = PinwheelSpecimenGeometry(cParams)
            obj.init(cParams);
            obj.compute();
        end

        function li = computeInnerLength(obj)
            [p,pIn,d] = obj.computeCutLine(1);
            [pPrev,pInPrev,dPrev] = obj.computeCutLine(obj.nSides);
            [~,pInNext,dNext] = obj.computeCutLine(2);
            pIn0 = obj.intersect(p,d,pInPrev,dPrev);
            pIn1 = obj.intersect(p,d,pInNext,dNext);
            li = norm(pIn1-pIn0);
            assert(isfinite(li) && all(isfinite([pPrev,pIn])));
        end

    end

    methods (Static, Access = public)

        function cP = hexagonParams(nCells)
            cP.nSides     = 3;
            cP.sideLength = 10;
            cP.t          = 0.18*cP.sideLength;
            cP.theta      = pi/33;
            cP.s          = 0.125;
            cP.hingeWidth = 0.25;
            cP.nCells     = nCells;
            cP.tiling     = 'hexagons';
        end

        function cP = hexagonPrintParams(nCells)
            cP.nSides     = 3;
            cP.sideLength = 40;
            cP.t          = 0.18*cP.sideLength;
            cP.theta      = pi/33;
            cP.s          = 0.8;
            cP.hingeWidth = 1.2;
            cP.nCells     = nCells;
            cP.tiling     = 'strip';
        end

        function cP = squareParams(nCells,tOverL)
            l = 20;
            t = tOverL*l;
            a = l/2;
            cP.nSides     = 4;
            cP.sideLength = l;
            cP.t          = t;
            cP.theta      = acos((a/2)/hypot(l/2,l/2-t)) - atan((l/2)/(l/2-t));
            cP.s          = 0.25;
            cP.hingeWidth = t;
            cP.nCells     = nCells;
            cP.tiling     = 'squares';
        end

    end

    methods (Access = protected)

        function w = computeCutWidth(obj)
            w = obj.s;
        end

        function shape = computeRefinementShape(obj)
            shape = obj.cuts;
        end

        function points = computeRefinePoints(obj)
            points = [obj.hingeTips;obj.latticeVertices];
            if obj.t <= 2*obj.hingeWidth
                points = [points;obj.cutRoots];
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.nSides     = cParams.nSides;
            obj.sideLength = cParams.sideLength;
            obj.t          = cParams.t;
            obj.theta      = cParams.theta;
            obj.s          = cParams.s;
            obj.hingeWidth = cParams.hingeWidth;
            obj.nCells     = cParams.nCells;
            obj.tiling     = cParams.tiling;
        end

        function compute(obj)
            [atomCuts,atomTips] = obj.computeAtomCuts();
            atomRoots = cell2mat(cellfun(@(c) c([1 4],:),atomCuts,'UniformOutput',false));
            [transforms,cells]  = obj.computeTiling();
            cuts = cell(numel(atomCuts)*numel(transforms),1);
            tips = cell(numel(transforms),1);
            roots    = cell(numel(transforms),1);
            vertices = cell(numel(transforms),1);
            ind  = 0;
            for iT = 1:numel(transforms)
                for iC = 1:numel(atomCuts)
                    ind = ind+1;
                    cuts{ind} = round(transforms{iT}(atomCuts{iC}),9);
                end
                tips{iT}     = transforms{iT}(atomTips);
                roots{iT}    = transforms{iT}(atomRoots);
                vertices{iT} = transforms{iT}(obj.computeAtomVertices());
            end
            domain = obj.mergePolygons(cells);
            voids  = obj.mergePolygons(obj.joinMirroredCuts(cuts));
            obj.cuts      = intersect(domain,voids);
            obj.shape     = subtract(domain,voids);
            obj.hingeTips = unique(round(vertcat(tips{:}),9),'rows');
            obj.cutRoots  = unique(round(vertcat(roots{:}),9),'rows');
            allVertices   = unique(round(vertcat(vertices{:}),9),'rows');
            [~,isOn]      = isinterior(obj.shape,allVertices(:,1),allVertices(:,2));
            obj.latticeVertices = allVertices(isOn,:);
            obj.center    = (min(allVertices)+max(allVertices))/2;
            obj.cornerLigament = obj.t;
        end

        function [cuts,tips] = computeAtomCuts(obj)
            n    = obj.nSides;
            cuts = cell(n,1);
            tips = zeros(2*n,2);
            for k = 1:n
                [p,pIn,d] = obj.computeCutLine(k);
                [~,pInNext,dNext] = obj.computeCutLine(mod(k,n)+1);
                normalNext = [-dNext(2),dNext(1)];
                tipLine  = pInNext + obj.hingeWidth*normalNext;
                tipOuter = obj.intersect(p,d,tipLine,dNext);
                tipInner = obj.intersect(pIn,d,tipLine,dNext);
                cuts{k}  = [p;tipOuter;tipInner;pIn];
                tips(2*k-1:2*k,:) = [tipOuter;tipInner];
            end
        end

        function [p,pIn,d] = computeCutLine(obj,k)
            n     = obj.nSides;
            V     = obj.computeAtomVertices();
            eOut  = V(mod(k,n)+1,:) - V(k,:);
            eIn   = V(k,:) - V(mod(k-2,n)+1,:);
            eOut  = eOut/norm(eOut);
            eIn   = eIn/norm(eIn);
            p     = V(k,:) - obj.t*eIn;
            pIn   = p - obj.s*eIn;
            c     = cos(obj.theta);
            sn    = sin(obj.theta);
            d     = [c*eOut(1)-sn*eOut(2), sn*eOut(1)+c*eOut(2)];
        end

        function V = computeAtomVertices(obj)
            L = obj.sideLength;
            switch obj.nSides
                case 3
                    V = [0 0; L 0; L/2 sqrt(3)*L/2];
                case 4
                    V = [0 0; L 0; L L; 0 L];
            end
        end

        function [transforms,cells] = computeTiling(obj)
            switch obj.tiling
                case 'hexagons'
                    [transforms,cells] = obj.computeHexagonTiling();
                case 'squares'
                    [transforms,cells] = obj.computeSquareTiling();
                case 'strip'
                    [transforms,cells] = obj.computeStripTiling();
            end
        end

        function atomMaps = computeHexagonAtomMaps(obj)
            atomMaps = cell(6,1);
            for j = 0:5
                if mod(j,2) == 0
                    atomMaps{j+1} = obj.rotation(j*pi/3);
                else
                    atomMaps{j+1} = obj.rotation((j+1)*pi/3)*diag([1,-1]);
                end
            end
        end

        function [transforms,cells] = computeStripTiling(obj)
            L = obj.sideLength;
            h = sqrt(3)*L/2;
            H = 2*obj.nCells*h;
            atomMaps = obj.computeHexagonAtomMaps();
            V = obj.computeAtomVertices();
            centroid = mean(V,1);
            transforms = {};
            cells      = {};
            for i = -2:ceil(H/(1.5*L))+2
                for jRow = -2:ceil(H/(2*h))+2
                    c = [1.5*L*i, 2*h*jRow + mod(i,2)*h];
                    for j = 1:6
                        p = centroid*atomMaps{j}' + c;
                        if p(1) > 0 && p(1) < H && p(2) > 0 && p(2) < H
                            map = @(x) x*atomMaps{j}' + c;
                            transforms{end+1,1} = map;
                            cells{end+1,1}      = map(V);
                        end
                    end
                end
            end
            obj.bottomGrips = [-L, H+L, 0];
            obj.topGrips    = [-L, H+L, H];
            obj.height      = H;
        end

        function [transforms,cells] = computeHexagonTiling(obj)
            L = obj.sideLength;
            h = sqrt(3)*L/2;
            n = obj.nCells;
            angles   = (0:5)'*pi/3;
            hexagon  = L*[cos(angles),sin(angles)];
            atomMaps = obj.computeHexagonAtomMaps();
            transforms = cell(6*n^2,1);
            cells      = cell(n^2,1);
            ind = 0;
            yBot = zeros(n,1);
            yTop = zeros(n,1);
            xCol = zeros(n,1);
            for i = 0:n-1
                xCol(i+1) = 1.5*L*i;
                yBot(i+1) = mod(i,2)*h - h;
                yTop(i+1) = yBot(i+1) + n*2*h;
                for jRow = 0:n-1
                    c = [xCol(i+1), yBot(i+1) + h + 2*h*jRow];
                    cells{i*n+jRow+1} = hexagon + c;
                    for j = 1:6
                        ind = ind+1;
                        transforms{ind} = @(x) x*atomMaps{j}' + c;
                    end
                end
            end
            obj.bottomGrips = [xCol-L/2, xCol+L/2, yBot];
            obj.topGrips    = [xCol-L/2, xCol+L/2, yTop];
            obj.height      = n*2*h;
        end

        function [transforms,cells] = computeSquareTiling(obj)
            L = obj.sideLength;
            n = obj.nCells;
            atomMaps = {diag([1,1]),diag([-1,1]),diag([1,-1]),diag([-1,-1])};
            unitCell = L*[-1 -1; 1 -1; 1 1; -1 1];
            transforms = cell(4*n^2,1);
            cells      = cell(n^2,1);
            ind = 0;
            for i = 0:n-1
                for jRow = 0:n-1
                    c = 2*L*[i,jRow];
                    cells{i*n+jRow+1} = unitCell + c;
                    for j = 1:4
                        ind = ind+1;
                        transforms{ind} = @(x) x*atomMaps{j}' + c;
                    end
                end
            end
            xMin = -L;
            xMax = 2*L*n - L;
            obj.bottomGrips = [xMin, xMax, -L];
            obj.topGrips    = [xMin, xMax, 2*L*n-L];
            obj.height      = 2*L*n;
        end

        function joined = joinMirroredCuts(obj,cuts)
            roots = cell2mat(cellfun(@(c) [c(1,:),c(4,:)],cuts,'UniformOutput',false));
            [~,~,group] = uniquetol(roots,1e-6,'ByRows',true,'DataScale',1);
            joined = cell(max(group),1);
            for i = 1:numel(joined)
                members = find(group == i);
                A = cuts{members(1)};
                if numel(members) == 2
                    B = cuts{members(2)};
                    joined{i} = [A(2,:);A(3,:);A(4,:);B(3,:);B(2,:);B(1,:)];
                else
                    d = (A(2,:)-A(1,:))/norm(A(2,:)-A(1,:));
                    outwards = 0.01*obj.sideLength*d;
                    A([1 4],:) = A([1 4],:) - outwards;
                    joined{i}  = A;
                end
            end
        end

        function merged = mergePolygons(obj,polygons)
            shapes = cell(numel(polygons),1);
            for i = 1:numel(polygons)
                shapes{i} = obj.toPolyshape(polygons{i});
            end
            merged = union([shapes{:}]);
        end

    end

    methods (Static, Access = private)

        function x = intersect(p,d,q,e)
            a = [d(:),-e(:)]\(q(:)-p(:));
            x = p + a(1)*d;
        end

        function R = rotation(a)
            R = [cos(a) -sin(a); sin(a) cos(a)];
        end

        function ps = toPolyshape(x)
            x = round(x,9);
            signedArea = sum(x(:,1).*x([2:end,1],2) - x([2:end,1],1).*x(:,2));
            if signedArea > 0
                x = flipud(x);
            end
            ps = polyshape(x(:,1),x(:,2),'KeepCollinearPoints',true);
        end

    end

end
