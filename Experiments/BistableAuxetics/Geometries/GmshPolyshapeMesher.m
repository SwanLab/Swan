classdef GmshPolyshapeMesher < handle

    properties (Access = private)
        shape
        refinementShape
        boundaryPoints
        refinePoints
        sizeMin
        sizeMax
        distMin
        distMax
        workDir
        name
        isBackground
    end

    methods (Access = public)

        function obj = GmshPolyshapeMesher(cParams)
            obj.init(cParams);
        end

        function mesh = compute(obj)
            geoFile = fullfile(obj.workDir,[obj.name,'.geo']);
            mshFile = fullfile(obj.workDir,[obj.name,'.msh']);
            obj.writeGeo(geoFile);
            obj.runGmsh(geoFile,mshFile);
            [s.coord,s.connec] = obj.readMsh2(mshFile);
            mesh = Mesh.create(s);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.shape           = cParams.shape;
            obj.refinementShape = cParams.refinementShape;
            obj.boundaryPoints  = cParams.boundaryPoints;
            obj.refinePoints    = cParams.refinePoints;
            obj.sizeMin         = cParams.sizeMin;
            obj.sizeMax         = cParams.sizeMax;
            obj.distMin         = cParams.distMin;
            obj.distMax         = cParams.distMax;
            obj.workDir         = cParams.workDir;
            obj.name            = cParams.name;
            obj.isBackground    = cParams.isBackground;
        end

        function writeGeo(obj,fileName)
            assert(obj.shape.NumRegions == 1,'The specimen must be connected');
            fid = fopen(fileName,'w');
            obj.writeHeader(fid);
            if obj.isBackground
                obj.writeBackground(fid);
            else
                obj.writeShape(fid);
            end
            fclose(fid);
        end

        function writeHeader(~,fid)
            fprintf(fid,'Mesh.MeshSizeFromPoints = 0;\n');
            fprintf(fid,'Mesh.MeshSizeExtendFromBoundary = 0;\n');
            fprintf(fid,'Mesh.MeshSizeFromCurvature = 0;\n');
            fprintf(fid,'Mesh.Algorithm = 6;\n');
        end

        function writeBackground(obj,fid)
            [xLim,yLim] = boundingbox(obj.shape);
            box = [xLim([1 2 2 1])',yLim([1 1 2 2])'];
            obj.writeLoop(fid,box,0,0,1);
            obj.writePlaneSurface(fid,1);
            ids = obj.writePoints(fid,obj.sampleBoundary(),size(box,1));
            obj.writeThreshold(fid,ids);
        end

        function ids = writePoints(~,fid,points,offset)
            ids = offset + (1:size(points,1));
            for i = 1:size(points,1)
                fprintf(fid,'Point(%d) = {%.12g,%.12g,0};\n',ids(i),points(i,1),points(i,2));
            end
        end

        function writePlaneSurface(~,fid,loops)
            loopStr = sprintf('%d,',loops);
            fprintf(fid,'Plane Surface(1) = {%s};\n',loopStr(1:end-1));
            fprintf(fid,'Physical Surface(1) = {1};\n');
        end

        function points = sampleBoundary(obj)
            points = zeros(0,2);
            for b = 1:numboundaries(obj.refinementShape)
                [x,y] = boundary(obj.refinementShape,b);
                for i = 1:numel(x)-1
                    points = [points;obj.sampleSegment([x(i),y(i)],[x(i+1),y(i+1)])];
                end
            end
        end

        function points = sampleSegment(obj,p,q)
            n = max(1,ceil(norm(q-p)/obj.sizeMin));
            s = (0:n-1)'/n;
            points = p + s.*(q-p);
        end

        function writeLoop(obj,fid,x,iP,iL,loop)
            n = size(x,1);
            obj.writePoints(fid,x,iP);
            for i = 1:n
                fprintf(fid,'Line(%d) = {%d,%d};\n',iL+i,iP+i,iP+mod(i,n)+1);
            end
            lines = sprintf('%d,',iL+(1:n));
            fprintf(fid,'Curve Loop(%d) = {%s};\n',loop,lines(1:end-1));
        end

        function writeShape(obj,fid)
            allPoints = zeros(0,2);
            for b = 1:numboundaries(obj.shape)
                [x,y] = boundary(obj.shape,b);
                [x,y] = obj.insertBoundaryPoints(x(1:end-1),y(1:end-1));
                n = size(allPoints,1);
                obj.writeLoop(fid,[x,y],n,n,b);
                allPoints = [allPoints;x,y];
            end
            obj.writePlaneSurface(fid,[find(~ishole(obj.shape));find(ishole(obj.shape))]);
            obj.writeSizeField(fid,allPoints);
        end

        function [x,y] = insertBoundaryPoints(obj,x,y)
            tol = 1e-7;
            for i = 1:size(obj.boundaryPoints,1)
                p = obj.boundaryPoints(i,:);
                if min(hypot(x-p(1),y-p(2))) > tol
                    next = [2:numel(x),1];
                    dLen = hypot(x(next)-x,y(next)-y);
                    gap  = hypot(x-p(1),y-p(2)) + hypot(x(next)-p(1),y(next)-p(2)) - dLen;
                    iSeg = find(gap < tol,1);
                    if ~isempty(iSeg)
                        x = [x(1:iSeg);p(1);x(iSeg+1:end)];
                        y = [y(1:iSeg);p(2);y(iSeg+1:end)];
                    end
                end
            end
        end

        function writeSizeField(obj,fid,points)
            [isTip,~] = ismembertol(points,obj.refinePoints,1e-6,'ByRows',true,'DataScale',1);
            assert(any(isTip),'No hinge tip found on the boundary');
            obj.writeThreshold(fid,find(isTip));
        end

        function writeThreshold(obj,fid,pointIds)
            ids = sprintf('%d,',pointIds);
            fprintf(fid,'Field[1] = Distance;\n');
            fprintf(fid,'Field[1].PointsList = {%s};\n',ids(1:end-1));
            fprintf(fid,'Field[2] = Threshold;\n');
            fprintf(fid,'Field[2].InField = 1;\n');
            fprintf(fid,'Field[2].SizeMin = %.12g;\n',obj.sizeMin);
            fprintf(fid,'Field[2].SizeMax = %.12g;\n',obj.sizeMax);
            fprintf(fid,'Field[2].DistMin = %.12g;\n',obj.distMin);
            fprintf(fid,'Field[2].DistMax = %.12g;\n',obj.distMax);
            fprintf(fid,'Background Field = 2;\n');
        end

        function runGmsh(~,geoFile,mshFile)
            cmd = sprintf('LD_LIBRARY_PATH= gmsh -2 -format msh2 -v 2 "%s" -o "%s"',geoFile,mshFile);
            [status,out] = system(cmd);
            assert(status == 0,'Gmsh failed: %s',out);
        end

        function [coord,connec] = readMsh2(~,fileName)
            lines  = splitlines(fileread(fileName));
            iNodes = find(strcmp(lines,'$Nodes'),1);
            nNodes = str2double(lines{iNodes+1});
            nodes  = sscanf(strjoin(lines(iNodes+2:iNodes+1+nNodes),' '),'%f',[4,nNodes])';
            iElems = find(strcmp(lines,'$Elements'),1);
            nElems = str2double(lines{iElems+1});
            connec = zeros(nElems,3);
            nTri   = 0;
            for i = 1:nElems
                row = sscanf(lines{iElems+1+i},'%d')';
                if row(2) == 2
                    nTri = nTri+1;
                    connec(nTri,:) = row(end-2:end);
                end
            end
            connec = connec(1:nTri,:);
            map = zeros(max(nodes(:,1)),1);
            map(nodes(:,1)) = 1:nNodes;
            connec = map(connec);
            [used,~,newIds] = unique(connec(:));
            connec = reshape(newIds,size(connec));
            coord  = nodes(used,2:3);
            a = coord(connec(:,2),:) - coord(connec(:,1),:);
            b = coord(connec(:,3),:) - coord(connec(:,1),:);
            isCW = a(:,1).*b(:,2) - a(:,2).*b(:,1) < 0;
            connec(isCW,[2 3]) = connec(isCW,[3 2]);
        end

    end

end
