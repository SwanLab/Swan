classdef MultilevelMesh < handle
    
    properties (Access = public)
        coarseMeshes
        interpolator
        mesh
    end
    
    properties (Access = private)
        nX
        nY
        nZ
        nLevel
    end
    
    methods (Access = public)
        
        function obj = MultilevelMesh(cParams)
            obj.init(cParams)         
            obj.createCoarseMeshesAndInterpolators();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)        
            obj.nX     = cParams.nX;
            obj.nY     = cParams.nY;
            obj.nZ     = cParams.nZ;
            obj.nLevel = cParams.nLevel;
        end
        
        function createCoarseMeshesAndInterpolators(obj)
           % nMesh   = obj.createCoarseMesh();
           length = 1;
           height = 1;
           width = 1;

           nMesh = obj.TetraMesh(length,height,width,obj.nX,obj.nY,obj.nZ);

           nInterp = obj.createInterpolator(nMesh);
           for iLevel = 1:obj.nLevel
              obj.coarseMeshes{iLevel} = nMesh;
              obj.interpolator{iLevel} = nInterp;               
              m                        = obj.coarseMeshes{iLevel};
              int                      = obj.interpolator{iLevel};
              nMesh                    = obj.remesh(m,int);
              nInterp                  = obj.createInterpolator(nMesh);            
           end
         obj.mesh = nMesh;
        end
        
       function meshFine = remesh(obj,mesh,I)
            meshCoarse = mesh;
            meshCoord  = I * meshCoarse.coord;
            meshConnec = delaunay(meshCoord);
            s.coord    = meshCoord;
            s.connec   = meshConnec;
            meshFine   = Mesh(s);
       end        
        
       function I = createInterpolator(obj,mesh)
            p = mesh.coord;
            t = mesh.connec;

            n = size(p,1);
            q = size(t,1);
            T = sparse(eye(n,n));
            tnew = []; j = 1;
            p_ori = p;
            for i = 1:q % this will add all the midpoints into p
                tcurr = t(i,:);
                pmid = [
                    (p(tcurr(1),:) + p(tcurr(2),:)) / 2;
                    (p(tcurr(2),:) + p(tcurr(3),:)) / 2;
                    (p(tcurr(3),:) + p(tcurr(1),:)) / 2;
                    (p(tcurr(1),:) + p(tcurr(4),:)) / 2;
                    (p(tcurr(2),:) + p(tcurr(4),:)) / 2;
                    (p(tcurr(3),:) + p(tcurr(4),:)) / 2;
                    ];
                p = [p; pmid];
            end

            [~,ia] = unique(p,'rows','stable');
            Ia = ia(n+1:end);
            ias = ia(n+1:end) - n ;
            potential_tri = ceil(ias./6);
            d = 1;
            midpt_curr = [];

            for i = potential_tri' % now need to loop thru ia and find the triangle that
                %corresponds to this midpoint
                tcurr = t(i,:);
                midpt_curr(1,:) = p(Ia(d),:);

                pmid = [
                    (p(tcurr(1),:) + p(tcurr(2),:)) / 2;
                    (p(tcurr(2),:) + p(tcurr(3),:)) / 2;
                    (p(tcurr(3),:) + p(tcurr(1),:)) / 2;
                    (p(tcurr(1),:) + p(tcurr(4),:)) / 2;
                    (p(tcurr(2),:) + p(tcurr(4),:)) / 2;
                    (p(tcurr(3),:) + p(tcurr(4),:)) / 2;
                    ];

                if midpt_curr(1,:) == pmid(1,:)
                    T(n + 1, [tcurr(1),tcurr(2)]) = 1/2;
                elseif midpt_curr(1,:) == pmid(2,:)
                    T(n + 1, [tcurr(2),tcurr(3)]) = 1/2;
                elseif midpt_curr(1,:) == pmid(3,:)
                    T(n + 1, [tcurr(1),tcurr(3)]) = 1/2;
                elseif midpt_curr(1,:) == pmid(4,:)
                    T(n + 1, [tcurr(1),tcurr(4)]) = 1/2;
                elseif midpt_curr(1,:) == pmid(5,:)
                    T(n + 1, [tcurr(2),tcurr(4)]) = 1/2;
                elseif midpt_curr(1,:) == pmid(6,:)
                    T(n + 1, [tcurr(3),tcurr(4)]) = 1/2;
                end
                n = n + 1;
                d = d + 1;
            end
            I = T;
        end
        
       
        function m = createCoarseMesh(obj)
        
            % Generate coordinates
            x1 = linspace(0,2,obj.nX);
            x2 = linspace(0,1,obj.nY);
            x3 = linspace(0,1,obj.nZ);
            % Create the grid
            [xv,yv,zv] = meshgrid(x1,x2,x3);
            % Triangulate the mesh to obtain coordinates and connectivities
            meshConnec = delaunayn([xv,yv,zv]);
            %[F,V] = mesh2tri(xv,yv,zeros(length(xv)),'x');

            s.coord = V(:,1:2);
            s.connec = F;
            m = Mesh(s);        
        end

        function m = createCoarseMesh3D(obj)
            filename   = 'test3d_tetrahedra';
            a.fileName = filename;
            femD       = FemDataContainer(a);
            m         = femD.mesh;
        end

        function m = TetraMesh(obj,length, height, width, nx, ny, nz)
           cMeshGlobal =  obj.HexaMesh(length, height, width, nx, ny, nz);
            Xiso     =  [-1 ,-1, -1;...
                1, -1, -1;...
                1, 1, -1;...
                -1, 1, -1;...
                -1, -1, 1;...
                1, -1, 1;...
                1, 1, 1;...
                -1, 1, 1;];
            connecIso    = delaunay(Xiso);
            ss.coord     = Xiso;
            ss.connec    = connecIso;
            localMesh    = Mesh(ss);
            nelem        = cMeshGlobal.nelem;
            bCutConnec   = cMeshGlobal.connec;
            connecIso    = localMesh.connec;
            nElemIso     = size(connecIso,1);
            nnodeSubMesh = size(connecIso,2);
            subConnec    = bCutConnec(:,connecIso');
            subConnec    = reshape(subConnec',[nnodeSubMesh,nelem*nElemIso])';
            sss.coord    = cMeshGlobal.coord;
            sss.connec   = subConnec;
            subMesh      = Mesh(sss);
            m = subMesh.computeCanonicalMesh();
        end

        function mesh = HexaMesh(obj,length, height, width, nx, ny, nz)
            % Coord
            x = linspace(0, length, nx+1);
            y = linspace(0, height, ny+1);
            z = linspace(0, width,  nz+1);
            [X,Y,Z] = meshgrid(x,y,z);
            npnod = size(X,1)*size(X,2)*size(X,3);
            Xr = reshape(X, npnod,1);
            Yr = reshape(Y, npnod,1);
            Zr = reshape(Z, npnod,1);
            coor = [Xr, Yr, Zr];
            
            % Connec
            glob = [];
            next_x = ny + 1;
            next_z = (nx+1)*(ny+1);
            for iZ = 1: nz
                z_toadd = (nx+1)*(ny+1) * (iZ-1);
                for iX = 1: nx
                    x_toadd = z_toadd + next_x*(iX-1);
                    for iY = 1: ny
                        y_line = x_toadd + [iY, iY + 1];
                        z0_plane = [y_line, flip(next_x + y_line) ];
                        nods = [z0_plane, next_z + z0_plane];
                        connec = [nods(2), nods(3), nods(7), nods(6), ...
                            nods(1), nods(4), nods(8), nods(5)];
                        glob = [glob; connec];
                    end
                end
            end
            
            s.coord = coor;
            s.connec = glob;
            
            mesh = Mesh(s);
        end
    end
end