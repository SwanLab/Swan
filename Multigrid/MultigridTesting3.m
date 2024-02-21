classdef MultigridTesting3 < handle
    
    
    properties (Access = private)
        nDimf
        boundaryConditionsFine
        fineMaterial
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
        FineK
        Lchol
        FineKred
        coarseMesh
        CoarseK
        CoarseKred
        coarseMaterial
        boundaryConditionsCoarse
        fineMeshCoord
        fineMeshConnec
        fineMesh
        FineDispFun
        CoarseDispFun
        FineFred
        CoarseFred
        data
        
        nMesh
        mesh
        I
        material
        boundaryConditions
        dispFun
        Kred
        Fred
        
    end
    
    methods (Access = public)
        
        function obj = MultigridTesting3()
            close all;
            addpath(genpath(fileparts(mfilename('fullpath'))))
            obj.init()
            obj.createCoarseMesh(1)
            obj.createBoundaryConditions(0)
            obj.createMaterial(0)
            obj.computeKred(0)
            obj.computeFred(0)
            
            for i = 1:obj.nMesh-1
                
                obj.createMatrixInterpolation(i)
                obj.createMesh(i)
                obj.createBoundaryConditions(i)
                obj.createMaterial(i)
                obj.computeKred(i)
                obj.computeFred(i)
                
            end
             obj.createData();

        end

        function r = getdata(obj)
            r = obj.data;
        end

        function r = getBC(obj)
            r = obj.boundaryConditions;
        end

        function r = getMesh(obj)
            r = obj.mesh;
        end
    end

    methods (Access = private)
        
        function init(obj)
            obj.nDimf = 2;
            obj.nbasis = 20;
            obj.functionType = 'P1';
            obj.nMesh = 2;
        end

        function createCoarseMesh(obj,i)
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
            obj.mesh{i} = Mesh(s);
        end

        function createMatrixInterpolation(obj,nMesh)
            p = obj.mesh{nMesh}.coord;
            t = obj.mesh{nMesh}.connec;

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
            obj.I{nMesh} = T;
        end
        
        function createMesh(obj,i)
            meshCoord = obj.I{i} * obj.mesh{i}.coord;
            meshConnec = delaunayn(meshCoord);
            s.coord = meshCoord;
            s.connec = meshConnec;
            obj.mesh{i+1} = Mesh(s);
        end
        
        function createBoundaryConditions(obj,i)
            rawBc    = obj.createRawBoundaryConditions(i);
            dim = getFunDims(obj,i);
            rawBc.ndimf = dim.ndimf;
            rawBc.ndofs = dim.ndofs;
            s.mesh  = obj.fineMesh;
            s.scale = 'MACRO';
            s.bc    = {rawBc};
            s.ndofs = dim.ndofs;
            bc = BoundaryConditions(s);
            bc.compute();
            obj.boundaryConditions{i+1} = bc;
        end
        
        function dim = getFunDims(obj,i)
            s.fValues = obj.mesh{i+1}.coord;
            s.mesh = obj.mesh{i+1};
            disp = P1Function(s);
            d.ndimf  = disp.ndimf;
            d.nnodes = size(disp.fValues, 1);
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = obj.mesh{i+1}.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end
        
        function bc = createRawBoundaryConditions(obj,i)
            dirichletNodes = abs(obj.mesh{i+1}.coord(:,1)-0) < 1e-12;
            rightSide  = max(obj.mesh{i+1}.coord(:,1));
            isInRight = abs(obj.mesh{i+1}.coord(:,1)-rightSide)< 1e-12;
            isInMiddleEdge = abs(obj.mesh{i+1}.coord(:,2)-1.5) < 0.1;
            %forceNodes = isInRight & isInMiddleEdge;
            forceNodes = isInRight;
            nodes = 1:obj.mesh{i+1}.nnodes;
            bcDir = [nodes(dirichletNodes)';nodes(dirichletNodes)'];
            nodesdir=size(nodes(dirichletNodes),2);
            bcDir(1:nodesdir,end+1) = 1;
            bcDir(nodesdir+1:end,end) = 2;
            bcDir(:,end+1)=0;
            bc.dirichlet = bcDir;
            bc.pointload(:,1) = nodes(forceNodes);
            bc.pointload(:,2) = 2;
            bc.pointload(:,3) = -1/length(forceNodes);
        end
        
        function createMaterial(obj,i)
            s.mesh = obj.mesh{i+1};
            s.type = 'ELASTIC';
            s.scale = 'MACRO';
            ngaus = 1;
            Id = ones(obj.mesh{i+1}.nelem,ngaus);
            s.ptype = 'ELASTIC';
            s.pdim  = '2D';
            s.nelem = obj.mesh{i+1}.nelem;
            s.mesh  = obj.mesh{i+1};
            s.kappa = .9107*Id;
            s.mu    = .3446*Id;
            mat = Material.create(s);
            mat.compute(s);
            obj.material{i+1} = mat;
        end
        
        function computeKred(obj,i)
            obj.dispFun{i+1} = P1Function.create(obj.mesh{i+1}, obj.nDimf);
            K    = obj.computeStiffnessMatrix(obj.mesh{i+1},obj.material{i+1},obj.dispFun{i+1});
            obj.Kred{i+1} = obj.boundaryConditions{i+1}.fullToReducedMatrix(K);
        end

        function k = computeStiffnessMatrix(obj,mesh,material,displacementFun)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = mesh;
            s.fun      = displacementFun;
            % s.test      = displacementFun;
            % s.trial      = displacementFun;
            s.material = material;
            lhs = LHSintegrator.create(s);
            k   = lhs.compute();
        end 
        
        function computeFred(obj,i)
            RHS  = obj.createRHS(obj.mesh{i+1},obj.dispFun{i+1},obj.boundaryConditions{i+1});
            Fext = RHS.compute();
            obj.Fred{i+1} = obj.boundaryConditions{i+1}.fullToReducedVector(Fext);
        end
        
        function RHS = createRHS(~,mesh,dispFun,boundaryConditions)
            dim.ndimf  = dispFun.ndimf;
            dim.nnodes = size(dispFun.fValues, 1);
            dim.ndofs  = dim.nnodes*dim.ndimf;
            dim.nnodeElem = mesh.nnodeElem; % should come from interp..
            dim.ndofsElem = dim.nnodeElem*dim.ndimf;
            c.dim=dim;
            c.mesh=mesh;
            c.BC = boundaryConditions;
            RHS    = RHSintegrator_ElasticMacro(c);
        end
        
        function createData(obj)
            
            for i = 1:obj.nMesh
                obj.data(i).p = obj.mesh{i}.coord;
                obj.data(i).t = obj.mesh{i}.connec;
                if i < obj.nMesh
                    obj.data(i).T = obj.I{i};
                    obj.data(i).R = obj.I{i}';
                end
                obj.data(i).A = obj.Kred{i};
                obj.data(i).b = obj.Fred{i};
            end
            
        end

    end

end