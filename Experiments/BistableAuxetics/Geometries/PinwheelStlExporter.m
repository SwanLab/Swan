classdef PinwheelStlExporter < handle

    properties (Access = private)
        shape
        thickness
        scale
        fileName
    end

    methods (Access = public)

        function obj = PinwheelStlExporter(cParams)
            obj.shape     = cParams.shape;
            obj.thickness = cParams.thickness;
            obj.scale     = cParams.scale;
            obj.fileName  = cParams.fileName;
        end

        function info = export(obj)
            [nodes,faces] = obj.computeSolid();
            stlwrite(triangulation(faces,nodes),obj.fileName);
            info = obj.computeChecks(nodes,faces);
        end

    end

    methods (Access = private)

        function [nodes,faces] = computeSolid(obj)
            tri   = obj.triangulateShape();
            p     = obj.scale*tri.Points;
            cap   = tri.ConnectivityList;
            nP    = size(p,1);
            edges = freeBoundary(tri);
            walls = [edges(:,1),edges(:,2),edges(:,2)+nP; ...
                     edges(:,1),edges(:,2)+nP,edges(:,1)+nP];
            nodes = [p,zeros(nP,1);p,obj.thickness*ones(nP,1)];
            faces = [cap(:,[1 3 2]);cap+nP;walls];
        end

        function tri = triangulateShape(obj)
            tri = triangulation(obj.shape);
            p   = tri.Points;
            f   = tri.ConnectivityList;
            cen = (p(f(:,1),:)+p(f(:,2),:)+p(f(:,3),:))/3;
            tri = triangulation(f(isinterior(obj.shape,cen(:,1),cen(:,2)),:),p);
        end

        function info = computeChecks(obj,nodes,faces)
            e = sort([faces(:,[1 2]);faces(:,[2 3]);faces(:,[3 1])],2);
            [~,~,ic] = unique(e,'rows');
            info.isClosed = all(accumarray(ic,1) == 2);
            a = nodes(faces(:,1),:);
            b = nodes(faces(:,2),:);
            c = nodes(faces(:,3),:);
            info.volume   = sum(sum(a.*cross(b,c,2),2))/6;
            info.expected = area(obj.shape)*obj.scale^2*obj.thickness;
            info.size     = max(nodes)-min(nodes);
            info.nFaces   = size(faces,1);
        end
    end

end
