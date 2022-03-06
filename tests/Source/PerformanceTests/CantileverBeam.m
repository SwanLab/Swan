classdef CantileverBeam < handle
    
    properties (Access = public)
        connec
        coords
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
            obj.coords = coords;
        end

        function connec = computeConnec2D(obj, xdiv, ydiv)
            connec = [];
            for j = 0:1:ydiv-1
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
            z = linspace(0, obj.height, ydiv+1);
            [X,Y,Z] = meshgrid(x,y,z);
            npnod = size(X,1)*size(X,2)*size(X,3);
            Xr = reshape(X, npnod,1);
            Yr = reshape(Y, npnod,1);
            Zr = reshape(Z, npnod,1);
            coords = [Xr, Yr, Zr];
            obj.coords = coords;
        end

        function connec = computeConnec3D(obj, xdiv, ydiv)
            connec = [];
            for z = 0:1:ydiv-1
                for j = 0:1:ydiv-1
                    for i = 1:1:xdiv
                        addZ = (xdiv+1)*(ydiv+1)*z;
                        addZ1 = (xdiv+1)*(ydiv+1)*(z+1);
                        node1 = j*(xdiv+1)+i + addZ;
                        node2 = j*(xdiv+1)+i+1 + addZ;
                        node3 = (j+1)*(xdiv+1)+i + addZ;
                        node4 = (j+1)*(xdiv+1)+i+1 + addZ;
                        node5 = j*(xdiv+1)+i + addZ1;
                        node6 = j*(xdiv+1)+i+1 + addZ1;
                        node7 = (j+1)*(xdiv+1)+i + addZ1;
                        node8 = (j+1)*(xdiv+1)+i+1 + addZ1;
                        elem = [node1, node2, node3, node4, ...
                                node5, node6, node7, node8];
                        connec = [connec; elem];
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