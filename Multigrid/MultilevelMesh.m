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
        Length
        Height
        Width
        ndim
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
            obj.Length = cParams.length;
            obj.Height = cParams.height;
            obj.Width  = cParams.width ;
            obj.ndim   = cParams.ndim;
        end
        
        function createCoarseMeshesAndInterpolators(obj)
          
           nMesh   = obj.createCoarseMesh();

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
            meshFine   = Mesh.create(s);
       end        
        
       function I = createInterpolator(obj,mesh)
           s.meshType = mesh.type;
           s.mesh     = mesh;

           inter = Interpolator.create(s);
           I     = inter.I;
        end
        
       
        function m = createCoarseMesh(obj)
            
            if obj.ndim == 2
            %  % %   Triangular Mesh
            %     x1 = linspace(0,obj.Length,obj.nX);
            %     x2 = linspace(0,obj.Height,obj.nY);
            %     [xv,yv] = meshgrid(x1,x2);
            %     [F,V] = mesh2tri(xv,yv,zeros(length(xv)),'x');
            %     s.coord = V(:,1:2);
            %     s.connec = F;
            % %    m = Mesh.create(s);
                m = TriangleMesh(obj.Length, obj.Height, obj.nX, obj.nY);

                % % Quadrilater Mesh
          %       m = QuadMesh(obj.Length, obj.Height, obj.nX, obj.nY);
            elseif obj.ndim == 3
                m = TetraMesh(obj.Length,obj.Height,obj.Width,obj.nX,obj.nY,obj.nZ);
               
            end
        end

        % function m = createCoarseMesh3D(obj)
        %     filename   = 'test3d_tetrahedra';
        %     a.fileName = filename;
        %     femD       = FemDataContainer(a);
        %     m         = femD.mesh;
        % end

%         function m = TetraMesh(obj,length, height, width, nx, ny, nz)
%            cMeshGlobal =  obj.HexaMesh(length, height, width, nx, ny, nz);
%             Xiso     =  [-1 ,-1, -1;...
%                 1, -1, -1;...
%                 1, 1, -1;...
%                 -1, 1, -1;...
%                 -1, -1, 1;...
%                 1, -1, 1;...
%                 1, 1, 1;...
%                 -1, 1, 1;];
%             connecIso    = delaunay(Xiso);
%             ss.coord     = Xiso;
%             ss.connec    = connecIso;
%             localMesh    = Mesh(ss);
%             nelem        = cMeshGlobal.nelem;
%             bCutConnec   = cMeshGlobal.connec;
%             connecIso    = localMesh.connec;
%             nElemIso     = size(connecIso,1);
%             nnodeSubMesh = size(connecIso,2);
%             subConnec    = bCutConnec(:,connecIso');
%             subConnec    = reshape(subConnec',[nnodeSubMesh,nelem*nElemIso])';
%             sss.coord    = cMeshGlobal.coord;
%             sss.connec   = subConnec;
%             subMesh      = Mesh(sss);
%             m = subMesh.computeCanonicalMesh();
%         end
% 
%         function mesh = HexaMesh(obj,length, height, width, nx, ny, nz)
%             % Coord
%             x = linspace(0, length, nx+1);
%             y = linspace(0, height, ny+1);
%             z = linspace(0, width,  nz+1);
%             [X,Y,Z] = meshgrid(x,y,z);
%             npnod = size(X,1)*size(X,2)*size(X,3);
%             Xr = reshape(X, npnod,1);
%             Yr = reshape(Y, npnod,1);
%             Zr = reshape(Z, npnod,1);
%             coor = [Xr, Yr, Zr];
%             
%             % Connec
%             glob = [];
%             next_x = ny + 1;
%             next_z = (nx+1)*(ny+1);
%             for iZ = 1: nz
%                 z_toadd = (nx+1)*(ny+1) * (iZ-1);
%                 for iX = 1: nx
%                     x_toadd = z_toadd + next_x*(iX-1);
%                     for iY = 1: ny
%                         y_line = x_toadd + [iY, iY + 1];
%                         z0_plane = [y_line, flip(next_x + y_line) ];
%                         nods = [z0_plane, next_z + z0_plane];
%                         connec = [nods(2), nods(3), nods(7), nods(6), ...
%                             nods(1), nods(4), nods(8), nods(5)];
%                         glob = [glob; connec];
%                     end
%                 end
%             end
%             
%             s.coord = coor;
%             s.connec = glob;
%             
%             mesh = Mesh(s);
%         end
    end
end