classdef MultigridTesting < handle
    
    
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
    end
    
    methods (Access = public)
        
        function obj = MultigridTesting()
            close all;
            addpath(genpath(fileparts(mfilename('fullpath'))))
            obj.init()
            obj.mesh = obj.createMesh();
            rawBc    = obj.createRawBoundaryConditions();
            obj.boundaryConditions = obj.createBoundaryConditions(rawBc);
            obj.material = obj.createMaterial();
            
            dispFun = P1Function.create(obj.mesh, obj.nDimf);
            
            obj.K    = obj.computeStiffnessMatrix(obj.mesh,obj.material,dispFun);
            obj.Kred = obj.boundaryConditions.fullToReducedMatrix(obj.K);
%             R = sprand((obj.Kred));
             obj.D = diag(diag(obj.Kred));
%             obj.L = tril(obj.Kred,-1);
%             obj.Lt= obj.L';
%             obj.Lchol=ichol(obj.Kred);
            
            RHS  = obj.createRHS(obj.mesh,dispFun,obj.boundaryConditions);
            Fext = RHS.compute();
            Fred = obj.boundaryConditions.fullToReducedVector(Fext);
            
            x=obj.Kred\Fred;
            
            [obj.basisFvalues,obj.basisVec,obj.eigenVec] = obj.computeBasis(obj.K);
            
%             modalFun  = ModalFunction.create(obj.mesh,obj.basisFvalues,obj.functionType);
% 
%             LHSK   = obj.createLHSstiffnessGlobal(obj.mesh,obj.quad,obj.material,modalFun);
%             obj.Kmodal = LHSK.compute();
            
            [xCG,residualCG]   = obj.conjugateGradient(obj.Kred,Fred,x);

%             GD.x        = xCG; 
%             GD.residual = residualCG; 
%             GD.err      = errCG;
%             GD.errA     = errACG;
            

            obj.coarseMesh = obj.createCoarseMesh();
            dispCoarseFun = P1Function.create(obj.coarseMesh, obj.nDimf);
            obj.KCoarse    = obj.computeStiffnessMatrix(obj.coarseMesh,obj.material,dispCoarseFun);
            obj.KredCoarse = obj.boundaryConditions.fullToReducedMatrix(obj.K);
            RHSCoarse  = obj.createRHS(obj.coarseMesh,dispCoarseFun,obj.boundaryConditions);
            FextCoarse = RHSCoarse.compute();
            FredCoarse = obj.boundaryConditions.fullToReducedVector(FextCoarse);
            xCoarse=obj.KCoarse\FredCoarse;

            [xPCG,residualPCG] = obj.preconditionedConjugateGradient(obj.Kred,Fred,x);

%             PGDbasis20.x        = xPCG; 
%             PGDbasis20.residual = residualPCG; 
%             PGDbasis20.err      = errPCG;
%             PGDbasis20.errA     = errAPCG;            
            
        end
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.nDimf = 2;
            obj.nbasis = 20;
            obj.functionType = 'P1';
%             for i = 1:obj.nbasis
%                 obj.functionType(i) = funTyp;
%             end
        end
                
        function bc = createRawBoundaryConditions(obj)
            dirichletNodes = abs(obj.mesh.coord(:,1)-0) < 1e-12;
            rightSide  = max(obj.mesh.coord(:,1));
            isInRight = abs(obj.mesh.coord(:,1)-rightSide)< 1e-12;
            isInMiddleEdge = abs(obj.mesh.coord(:,2)-1.5) < 0.1;
            forceNodes = isInRight & isInMiddleEdge;
            nodes = 1:obj.mesh.nnodes;
            bcDir = [nodes(dirichletNodes)';nodes(dirichletNodes)'];
            nodesdir=size(nodes(dirichletNodes),2);
            bcDir(1:nodesdir,end+1) = 1;
            bcDir(nodesdir+1:end,end) = 2;
            bcDir(:,end+1)=0;
            bc.dirichlet = bcDir;
            bc.pointload(:,1) = nodes(forceNodes);
            bc.pointload(:,2) = 2;
            bc.pointload(:,3) = -1;
        end
        
        function bc = createBoundaryConditions(mesh,bcV)
            dim = getFunDims(mesh);
            bcV.ndimf = dim.ndimf;
            bcV.ndofs = dim.ndofs;
            s.mesh  = mesh;
            s.scale = 'MACRO';
            s.bc    = {bcV};
            s.ndofs = dim.ndofs;
            bc = BoundaryConditions(s);
            bc.compute();
        end
        
        function dim = getFunDims(obj)
            s.fValues = obj.mesh.coord;
            s.mesh = obj.mesh;
            disp = P1Function(s);
            d.ndimf  = disp.ndimf;
            d.nnodes = size(disp.fValues, 1);
            d.ndofs  = d.nnodes*d.ndimf;
            d.nnodeElem = obj.mesh.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end
        
        function material = createMaterial(obj)
            s.mesh = obj.mesh;
            s.type = 'ELASTIC';
            s.scale = 'MACRO';
            ngaus = 1;
            I = ones(obj.mesh.nelem,ngaus);
            s.ptype = 'ELASTIC';
            s.pdim  = '2D';
            s.nelem = obj.mesh.nelem;
            s.mesh  = obj.mesh;
            s.kappa = .9107*I;
            s.mu    = .3446*I;
            mat = Material.create(s);
            mat.compute(s);
            material = mat;
        end
        
        function [basis,basisVec,eigenVec] = computeBasis(obj,K)
            [eigenVec,D]=eigs(K,obj.nbasis,'smallestabs');
            psi = K*eigenVec;

            for i = 1:size(eigenVec,2)
                b = eigenVec(:,i);
                %b1 = obj.boundaryConditions.reducedToFullVector(b);
                basis{i} = reshape(b,2,[])';
                basisVec{i}= b;
                %a = psi(:,i);
                %a1=obj.boundaryConditions.reducedToFullVector(a);
                %psiD{i}=reshape(a1,2,[])';
            end

        end
        
        function [x,residual] = preconditionedConjugateGradient(obj,A,B,xsol)
            tol = 1e-6;
            n = length(B);
            x = zeros(n,1);
            r = B - A * x;
            %             z = ModalTesting.applyPreconditioner(M,r);
            % z = obj.applyPrecoditionerMultigrid(r);
            %             z=r-z;
            p = r;
            rzold = r' * r;
            iter = 0;

            hasNotConverged = true;

            while hasNotConverged
                
                Ap = A * p;
                alpha = rzold / (p' * Ap);
                x = x + alpha * p;
                r = r - alpha * Ap;
                if iter == 0 || hasPartiallyConverged == true 
                    ri = r;
                end
                %                 z = ModalTesting.applyPreconditioner(M,r)
                rznew = r' * r;

                %hasNotConverged = sqrt(rsnew) > tol;
                hasNotConverged = norm(r) > tol;
                hasPartiallyConverged = norm(r)/norm(ri) < 0.5;

                if hasPartiallyConverged
                    z = obj.applyPrecoditionerMultigrid(r,KredCoarse,FredCoarse);
                    r = z;
                end

                p = (rznew / rzold) * p;
                p = r + p;
                rzold = rznew;
                iter = iter + 1;
                residual(iter) = norm(r); %Ax - b
%                 err(iter)=norm(x-xsol);
%                 errAnorm(iter)=((x-xsol)')*A*(x-xsol);
            end
        end
        %
%                 function z = applyPreconditioner(obj,r)
%                     lhs=obj.Kmodal;
%                     phi=obj.eigenVec;
%                     r1=phi'*r;
%                     zP=lhs\r1;
%                     z=phi*zP;
%                     %z = r;
%                     z = r-z;
%                 end
%                 function z = applyPreconditioner(obj,r)
%                     z=obj.D\r;
%                 end
                 function z = applyPrecoditionerMultigrid(~,u,A,B)
                    %R = generateR(u);
                    j = 1;
                    for i=1:2:length(u)-1
                       R(j,i) = 1;
                       R(j,i+1) = 2;
                       R(j,i+2) = 1;
                       j = j+1;
                    end
                    R(:,i) = [];
                    z = 1/4 .* R * u;

                    r = z;

                    n = length(B);
                    x = zeros(n,1);
                    p = r;
                    rzold = r' * r;
                    iter = 0;
        
                    hasNotPartyallyConverged = true;
        
                    while hasNotPartyallyConverged
                        
                        Ap = A * p;
                        alpha = rzold / (p' * Ap);
                        x = x + alpha * p;
                        r = r - alpha * Ap;

                        rznew = r' * r;
        
                        %hasNotPartyallyConverged = ;
        
                        p = (rznew / rzold) * p;
                        p = r + p;
                        rzold = rznew;
                        iter = iter + 1;
                        residual(iter) = norm(r); %Ax - b
                    end

                    j = 1;
                    for i=1:2:length(u)
                       I(i,j) = 1;
                       I(i+1,j) = 2;
                       I(i+2,j) = 1;
                       j = j+1;
                    end   
                    z = 1/2 .* I .* u;                        
                 end

%                 function R = generateR(~,u)
%                     j = 1;
%                    for i=1:2:length(u)
%                        R(i,j) = 1;
%                        R(i+1,j) = 2;
%                        R(1+2,j) = 1;
%                        j = j+1;
%                    end        
%                 end
                 
    end
    
    methods (Access = public, Static)
        
        function mesh = createMesh()
            % Generate coordinates
            x1 = linspace(0,2,30);
            x2 = linspace(1,2,30);
            % Create the grid
            [xv,yv] = meshgrid(x1,x2);
            % Triangulate the mesh to obtain coordinates and connectivities
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');

            s.coord = V(:,1:2);
            s.connec = F;
            mesh = Mesh(s);
        end

        function coarsemesh = createCoarseMesh()
            % Generate coordinates
            x1 = linspace(0,2,15);
            x2 = linspace(1,2,15);
            % Create the grid
            [xv,yv] = meshgrid(x1,x2);
            % Triangulate the mesh to obtain coordinates and connectivities
            [F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');

            s.coord = V(:,1:2);
            s.connec = F;
            coarsemesh = Mesh(s);
        end
        
        function k = computeStiffnessMatrix(mesh,material,displacementFun)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = mesh;
            s.fun      = displacementFun;
            % s.test      = displacementFun;
            % s.trial      = displacementFun;
            s.material = material;
            lhs = LHSintegrator.create(s);
            k   = lhs.compute();
        end 
        
        function RHS = createRHS(mesh,dispFun,boundaryConditions)
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
        
        function LHSglobal = createLHSstiffnessGlobal(mesh,quad,material,modalFun)
            sL.material = material;
            sL.test = modalFun;
            sL.trial = modalFun;
            sL.mesh = mesh;
            sL.quadratureOrder = quad.order;
            LHSglobal = LHS_integratorStiffnessGlobal(sL);
        end
        
        function [x,residual] = conjugateGradient(LHS,RHS,xsol)
            tol = 1e-6;
            n = length(RHS);
            x = zeros(n,1);
            r = RHS - LHS * x;
            p = r;
            rsold = r' * r;
            iter = 0;

            hasNotConverged = true;

            while hasNotConverged
                Ap = LHS * p;
                alpha = rsold / (p' * Ap);
                x = x + alpha * p;
                r = r - alpha * Ap;
                rsnew = r' * r;

                %hasNotConverged = sqrt(rsnew) > tol;
                hasNotConverged = norm(LHS*x - RHS) > tol;

                p = r + (rsnew / rsold) * p;
                rsold = rsnew;
                iter = iter + 1;
                residual(iter) = norm(LHS*x - RHS); %Ax - b
%                 err(iter)=norm(x-xsol);
%                 errAnorm(iter)=((x-xsol)')*LHS*(x-xsol);
%                 f(iter)= 0.5*(x'*LHS*x)-RHS'*x;
            end
        end
        
    end
    
end

