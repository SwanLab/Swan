classdef MultilevelMesh < handle
    
    properties (Access = public)
        coarseMeshes
        interpolator
        mesh
    end
    
    properties (Access = private)
        nX
        nY
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
            obj.nLevel = cParams.nLevel;
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
            meshConnec = delaunayn(meshCoord);
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
                pmid = [(p(tcurr(1),:) + p(tcurr(2),:)) / 2;
                    (p(tcurr(2),:) + p(tcurr(3),:)) / 2;
                    (p(tcurr(3),:) + p(tcurr(1),:)) / 2];
                p = [p; pmid];
            end

            [~,ia] = unique(p,'rows','stable');
            Ia = ia(n+1:end);
            ias = ia(n+1:end) - n ;
            potential_tri = ceil(ias./3);
            d = 1;
            midpt_curr = [];

            for i = potential_tri' % now need to loop thru ia and find the triangle that
                %corresponds to this midpoint
                tcurr = t(i,:);
                midpt_curr(1,:) = p(Ia(d),:);

                pmid = [(p(tcurr(1),:) + p(tcurr(2),:)) / 2;
                    (p(tcurr(2),:) + p(tcurr(3),:)) / 2;
                    (p(tcurr(3),:) + p(tcurr(1),:)) / 2];

                if midpt_curr(1,:) == pmid(1,:)
                    T(n + 1, [tcurr(1),tcurr(2)]) = 1/2;
                elseif midpt_curr(1,:) == pmid(2,:)
                    T(n + 1, [tcurr(2),tcurr(3)]) = 1/2;
                elseif midpt_curr(1,:) == pmid(3,:)
                    T(n + 1, [tcurr(1),tcurr(3)]) = 1/2;
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
            % Create the grid
            [xv,yv] = meshgrid(x1,x2);
            % Triangulate the mesh to obtain coordinates and connectivities
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');

            s.coord = V(:,1:2);
            s.connec = F;
            m = Mesh(s);        
        end
        
    end
    
end