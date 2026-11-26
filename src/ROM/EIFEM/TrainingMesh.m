classdef TrainingMesh < handle

    properties (Access =private)
        meshSource
        ndim
        nelem
        bounds
        meshFile
    end

    methods (Access = public)

        function obj = TrainingMesh(s)
            obj.init(s);
        end

        function mS = create(obj)
            switch obj.meshSource
                case 'Imported'
                    mS = obj.TrainingMeshImported();
                case 'Generated'
                    switch obj.ndim
                        case 2
                            mS = obj.TrainingMesh2D();
                        case 3
                            mS = obj.TrainingMesh3D();
                        otherwise
                            error('Dimension must be 2 or 3.');
                    end
                otherwise
                    error('Unknown mesh source.');
            end
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = private)

        function init(obj,cParams)
            obj.meshSource = cParams.meshSource;
            obj.ndim       = cParams.ndim;
            obj.nelem      = cParams.nelem;
            obj.bounds     = cParams.bounds;
            obj.meshFile   = cParams.meshFile;
        end

        function m= TrainingMeshImported(obj)
            data = load(obj.meshFile);
            if isfield(data, 'EIFEoper')
                s.coord    = data.EIFEoper.MESH.COOR;
                s.connec   = data.EIFEoper.MESH.CN;
            elseif isfield(data, 'mesh')
                s.coord  = data.mesh.coord;
                s.connec = data.mesh.connec;
            else
                error('The file must contain EIFEoper or mesh.');
            end

            switch obj.ndim
                case 2
                    m = obj.moveEdges2D(s);
                case 3
                    m = obj.moveEdges3D(s);
                otherwise
                    error('Dimension must be 2 or 3.');
            end
        end


        function m = TrainingMesh2D(obj)
            n        = obj.nelem;
            x1       = linspace(0,obj.bounds(1),n(1));
            x2       = linspace(0,obj.bounds(2),n(2));
            [xv,yv]  = meshgrid(x1,x2);
            [F,V]    = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            m = obj.moveEdges2D(s);
        end


        function m = TrainingMesh3D(obj)
            n = obj.nelem;
            m = TetraMesh(obj.bounds(1),obj.bounds(2),obj.bounds(3),n(1),n(2),n(3));
            s.coord  = m.coord;
            s.connec = m.connec;
            m = obj.moveEdges3D(s);
        end

    end


    methods (Access = private, Static)

        function mS = moveEdges2D(s)
            maxC= max(s.coord);
            minC = min(s.coord);
            xmin = minC(1);
            xmax = maxC(1);
            ymin = minC(2);
            ymax = maxC(2);

            delta = 1e-9;
            s.coord(s.coord(:,1)== xmax & s.coord(:,2)==ymax,:) =...
                s.coord(s.coord(:,1)== xmax & s.coord(:,2)==ymax,:)+[-delta,-delta];
            s.coord(s.coord(:,1)== xmax & s.coord(:,2)==ymin,:) =...
                s.coord(s.coord(:,1)== xmax & s.coord(:,2)==ymin,:)+[-delta,+delta];
            s.coord(s.coord(:,1)== xmin & s.coord(:,2)==ymax,:) =...
                s.coord(s.coord(:,1)== xmin & s.coord(:,2)==ymax,:)+[+delta,-delta];
            s.coord(s.coord(:,1)== xmin & s.coord(:,2)==ymin,:) =...
                s.coord(s.coord(:,1)== xmin & s.coord(:,2)==ymin,:)+[+delta,+delta];
            mS = Mesh.create(s);
        end


        function mS = moveEdges3D(s)
            maxC= max(s.coord);
            minC = min(s.coord);

            delta=1E-9;
            s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:)-[delta,0,0];

            s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:)-[delta,0,0];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:)+[delta,0,0];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:)+[delta,0,0];

            s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:)-[0,0,delta];

            s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:)+[0,0,delta];

            s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:) =...
                s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:)-[0,0,delta];

            s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:) =...
                s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:)+[0,0,delta];

            s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==maxC(2),:) =...
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==maxC(2),:)-[0,delta,0];

            s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==minC(2),:) =...
                s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==minC(2),:)+[0,delta,0];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==maxC(2),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==maxC(2),:)-[0,delta,0];

            s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==minC(2),:) =...
                s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==minC(2),:)+[0,delta,0];

            mS = Mesh.create(s);
        end

    end

end