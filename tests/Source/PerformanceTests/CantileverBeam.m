classdef CantileverBeam < handle
    
    properties (Access = public)
        connec
    end
    
    properties (Access = private)
        dim
        len
        height
    end
    
    methods (Access = public)
        
        function obj = CantileverBeam(cParams)
            obj.init(cParams)
        end

        function mesh = create(obj, xdiv, ydiv)
            switch obj.dim
                case '2D'
                    mesh = obj.create2Dcantilever(xdiv, ydiv);
                case '3D'
                    mesh = obj.create3Dcantilever(xdiv, ydiv);
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.dim    = cParams.dim;
            obj.len    = cParams.len;
            obj.height = cParams.height;
        end

        function mesh = create2Dcantilever(obj, xdiv, ydiv)
            coords = obj.computeCoords2D(xdiv, ydiv);
            connec = obj.computeConnec2D(xdiv, ydiv);
            mesh   = obj.createMesh(coords, connec);
        end

        function mesh = create3Dcantilever(obj, xdiv, ydiv)
            coords = obj.computeCoords3D(xdiv, ydiv);
            connec = obj.computeConnec3D(xdiv, ydiv);
            mesh   = obj.createMesh(coords, connec);
        end

        function coords = computeCoords2D(obj, xdiv, ydiv)
            x = linspace(0, obj.len,    xdiv+1);
            y = linspace(0, obj.height, ydiv+1);
            [X,Y,Z] = meshgrid(x,y,0);
            fvc = surf2patch(X,Y,Z,'triangles');
            fvc.vertices(:,3) = []; % 2D
            coords = fvc.vertices;
        end

        function connec = computeConnec2D(obj, xdiv, ydiv)
            connec = [];
            for j = 0:1:ydiv
                for i = 1:1:xdiv
                    node1 = j*(xdiv+1)+i;
                    node2 = j*(xdiv+1)+i+1;
                    node3 = (j+1)*(xdiv+1)+i;
                    node4 = (j+1)*(xdiv+1)+i+1;
                    elem = [node1, node2, node3, node4];
                    connec = [connec;elem];
                end
            end
            obj.connec = connec;
        end

        function coords = computeCoords3D(obj, xdiv, ydiv)
            x = linspace(0, obj.len,    xdiv+1);
            y = linspace(0, obj.height, ydiv+1);
            [X,Y,Z] = meshgrid(x,y,y);
            fvc = surf2patch(X,Y,Z,'triangles');
            fvc.vertices(:,3) = []; % 2D
            coords = fvc.vertices;
        end

        function connec = computeConnec3D(obj, xdiv, ydiv)
            connec = [];
            for z = 1:1:ydiv
                for j = 0:1:ydiv
                    for i = 1:1:xdiv
                        node1 = j*(xdiv+1)+i; % should include z
                        node2 = j*(xdiv+1)+i+1;
                        node3 = (j+1)*(xdiv+1)+i;
                        node4 = (j+1)*(xdiv+1)+i+1;
                        elem = [node1, node2, node3, node4];
                        connec(:,:,z) = [connec(:,:,z);elem];
                    end
                end
            end
            obj.connec = connec;
        end

        function mesh = createMesh(obj, coords, connec)
            m.coord = coords;
            m.connec = connec;
            mesh = Mesh(m);
        end
        
    end
    
end