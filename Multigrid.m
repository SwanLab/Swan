classdef Multigrid < handle
    
    
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
        Kred
        Lchol
    end
    
    methods (Access = public)
        
        function obj = Multigrid()
            close all;
            obj.init()
            obj.mesh = obj.createMesh();
            obj.boundaryConditions = obj.createBoundaryConditions();
            obj.material = obj.createMaterial();
            
            dispFun = P1Function.create(obj.mesh, obj.nDimf);
            
            K    = obj.computeStiffnessMatrix(obj.mesh,obj.material,dispFun);
%             obj.Kred = obj.boundaryConditions.fullToReducedMatrix(K);
%             R = sprand((obj.Kred));
%             obj.D = diag(diag(obj.Kred));
%             obj.L = tril(obj.Kred,-1);
%             obj.Lt= obj.L';
%             obj.Lchol=ichol(obj.Kred);
            
            RHS  = obj.createRHS(obj.mesh,dispFun,obj.boundaryConditions);
            %Fext = RHS.compute();
            %Fred = obj.boundaryConditions.fullToReducedVector(Fext);
            
            [obj.basisFvalues,obj.basisVec,obj.eigenVec] = obj.computeBasis(obj.Kred);
            
            modalFun  = ModalFunction.create(obj.mesh,obj.basisFvalues,obj.functionType);

            LHSK   = obj.createLHSstiffnessGlobal(obj.mesh,obj.quad,obj.material,modalFun);
            obj.Kmodal = LHSK.compute();
            
            [xCG,residualCG,errCG,errACG]   = obj.conjugateGradient(Kred,Fext,x);

            GD.x        = xCG; 
            GD.residual = residualCG; 
            GD.err      = errCG;
            GD.errA     = errACG;
            
            [xPCG,residualPCG,errPCG,errAPCG] = obj.preconditionedConjugateGradient(obj.Kred,Fred,x);

            PGDbasis20.x        = xPCG; 
            PGDbasis20.residual = residualPCG; 
            PGDbasis20.err      = errPCG;
            PGDbasis20.errA     = errAPCG;            
            
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
                
        function bc = createBoundaryConditions(obj)
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
                b1 = obj.boundaryConditions.reducedToFullVector(b);
                basis{i} = reshape(b1,2,[])';
                basisVec{i}= b1;
                a = psi(:,i);
                a1=obj.boundaryConditions.reducedToFullVector(a);
                %psiD{i}=reshape(a1,2,[])';
            end

        end
        
        function [x,residual,err,errAnorm] = preconditionedConjugateGradient(obj,A,B,xsol)
            tol = 1e-6;
            n = length(B);
            x = zeros(n,1);
            r = B - A * x;
            %             z = ModalTesting.applyPreconditioner(M,r);
            z = obj.applyPreconditioner(r);
            %             z=r-z;
            p = z;
            rzold = r' * z;
            iter = 0;

            hasNotConverged = true;

            while hasNotConverged
                Ap = A * p;
                alpha = rzold / (p' * Ap);
                x = x + alpha * p;
                r = r - alpha * Ap;
                %                 z = ModalTesting.applyPreconditioner(M,r);
                z = obj.applyPreconditioner(r);
                rznew = r' * z;

                %hasNotConverged = sqrt(rsnew) > tol;
                hasNotConverged = norm(r) > tol;

                p = z + (rznew / rzold) * p;
                rzold = rznew;
                iter = iter + 1;
                residual(iter) = norm(r); %Ax - b
                err(iter)=norm(x-xsol);
                errAnorm(iter)=((x-xsol)')*A*(x-xsol);
            end
        end
        %
                function z = applyPreconditioner(obj,r)
                    lhs=obj.Kmodal;
                    phi=obj.eigenVec;
                    r1=phi'*r;
                    zP=lhs\r1;
                    z=phi*zP;
                    %z = r;
                    z = r-z;
                end
        %
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
        
        function [x,residual,err,errAnorm] = conjugateGradient(LHS,RHS,xsol)
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

                p = r + 0*(rsnew / rsold) * p;
                rsold = rsnew;
                iter = iter + 1;
                residual(iter) = norm(LHS*x - RHS); %Ax - b
                 err(iter)=norm(x-xsol);
                errAnorm(iter)=((x-xsol)')*LHS*(x-xsol);
                f(iter)= 0.5*(x'*LHS*x)-RHS'*x;
            end
        end
        
    end
    
end

