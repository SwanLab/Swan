classdef MultigridTesting2 < handle
    
    
    properties (Access = public)
        nDimf
        mesh
        boundaryConditions
        material
        quad
        basisFvalues
        basisVec
        eigenVec
        functionType
        nbasis
        Kmodal
        Mmodal
        D
        L
        Lt
        K
        Lchol
        Kred
        coarseMesh
        KCoarse
        KredCoarse
        coarseMaterial
        boundaryConditionsCoarse
        I
        fineMeshCoord
        fineMeshConnec
        fineMesh
    end
    
    methods (Access = public)
        
        function obj = MultigridTesting2()
            close all;
            addpath(genpath(fileparts(mfilename('fullpath'))))
            obj.init()
            obj.createCoarseMesh()
            obj.createMatrixInterpolation()
            obj.createFineMesh()
            
        end
    end

    methods (Access = private)
        
        function init(obj)
            obj.nDimf = 2;
            obj.nbasis = 20;
            obj.functionType = 'P1';
        end

        function createCoarseMesh(obj)
            numero1 = 10;
            numero2 = 10;
            % Generate coordinates
            x1 = linspace(0,2,numero1);
            x2 = linspace(0,1,numero2);
            % Create the grid
            [xv,yv] = meshgrid(x1,x2);
            % Triangulate the mesh to obtain coordinates and connectivities
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');

            s.coord = V(:,1:2);
            s.connec = F;
            obj.coarseMesh = Mesh(s);
        end

        function createMatrixInterpolation(obj)
            p = obj.coarseMesh.coord;
            t = obj.coarseMesh.connec;

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
                    obj.I(n + 1, [tcurr(1),tcurr(2)]) = 1/2;
                elseif midpt_curr(1,:) == pmid(2,:)
                    obj.I(n + 1, [tcurr(2),tcurr(3)]) = 1/2;
                elseif midpt_curr(1,:) == pmid(3,:)
                    obj.I(n + 1, [tcurr(1),tcurr(3)]) = 1/2;
                end
                n = n + 1;
                d = d + 1;
            end
        end

        function createFineMesh(obj)
            obj.fineMeshCoord = obj.I * obj.coarseMesh.coord;
            obj.fineMeshConnec = delaunayn(obj.fineMeshCoord);
            s.coord = obj.fineMeshCoord;
            s.connec = obj.fineMeshConnec;
            obj.fineMesh = Mesh(s);
        end
            
    end

end