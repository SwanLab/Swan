classdef ModalTesting < handle

    properties (Access = public)

    end

    properties (Access = private)
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

    properties (Access = private)

    end

    methods (Access = public)

        function obj = ModalTesting()
            close all;
            obj.init()
            obj.mesh = obj.createMesh();
            rawBc    = obj.createRawBoundaryConditions();
            obj.boundaryConditions = obj.createBoundaryConditions(rawBc);
            obj.quad = Quadrature.set(obj.mesh.type);
            obj.quad.computeQuadrature('QUADRATIC');
            obj.material = obj.createMaterial();

            dispFun = P1Function.create(obj.mesh, obj.nDimf);

            K    = obj.computeStiffnessMatrix(obj.mesh,obj.material,dispFun);
            obj.Kred = obj.boundaryConditions.fullToReducedMatrix(K);
            R = sprand((obj.Kred));
%            obj.Kred=0.5*(R'*R);
            obj.D = diag(diag(obj.Kred));
            obj.L = tril(obj.Kred,-1);
            obj.Lt= obj.L';
            obj.Lchol=ichol(obj.Kred);

            M    = obj.computeMassMatrix(obj.mesh,dispFun);
            Mred = obj.boundaryConditions.fullToReducedMatrix(M);

            RHS  = obj.createRHS(obj.mesh,dispFun,obj.boundaryConditions);
            Fext = RHS.compute();
            Fred = obj.boundaryConditions.fullToReducedVector(Fext);

            x=obj.Kred\Fred;


            [obj.basisFvalues,obj.basisVec,obj.eigenVec] = obj.computeBasis(obj.Kred);

            modalFun  = ModalFunction.create(obj.mesh,obj.basisFvalues,obj.functionType);

            LHSK   = obj.createLHSstiffnessGlobal(obj.mesh,obj.quad,obj.material,modalFun);
            obj.Kmodal = LHSK.compute();
            %             LHSM   = obj.createLHSmassGlobal(obj.mesh,obj.quad,obj.material,modalFun);
            %             obj.Mmodal = LHSM.compute();

            %             precond = obj.createModalPreconditioner();

            [xCG,residualCG,errCG,errACG]   = obj.conjugateGradient(obj.Kred,Fred,x);

            GD.x        = xCG; 
            GD.residual = residualCG; 
            GD.err      = errCG;
            GD.errA     = errACG;

            [xPCG,residualPCG,errPCG,errAPCG] = obj.preconditionedConjugateGradient(obj.Kred,Fred,x);

            PGDbasis20.x        = xPCG; 
            PGDbasis20.residual = residualPCG; 
            PGDbasis20.err      = errPCG;
            PGDbasis20.errA     = errAPCG;

            plot(residualCG); hold on; plot(residualPCG);
            set(gca, 'YScale', 'log')
            legend('GD','PCG (ssor)')
            xlabel('iteration')
            ylabel('residual')

            figure
            plot(errCG); hold on; plot(errPCG)
            set(gca, 'YScale', 'log')
            legend('GD','PCG (ssor)')
            xlabel('iteration')
            ylabel('||x-x*||')

            figure
            plot(errACG); hold on; plot(errAPCG)
            set(gca, 'YScale', 'log')
            legend('GD','PCG (ssor)')
            xlabel('iteration')
            ylabel('||x-x*||_{A}')
        end

    end

    methods (Access = private)

        function init(obj)
            obj.nDimf   = 2;
            obj.nbasis = 10;
            funTyp='P1';
            for i=1:obj.nbasis
                obj.functionType{i}=funTyp;
            end

        end


        function bc = createRawBoundaryConditions(obj)
            dirichletNodes = abs(obj.mesh.coord(:,1)-0) < 1e-12;
            rightSide  = max(obj.mesh.coord(:,1));
            isInRight = abs(obj.mesh.coord(:,1)-rightSide)< 1e-12;
            isInMiddleEdge = abs(obj.mesh.coord(:,2)-0.5) < 0.1;
            forceNodes = isInRight & isInMiddleEdge;
            forceNodes = isInRight;
            nodes = 1:obj.mesh.nnodes;
            bcDir = [nodes(dirichletNodes)';nodes(dirichletNodes)'];
            % bcDir = [nodes(dirichletNodes)'];

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
            I = ones(obj.mesh.nelem,obj.quad.ngaus);
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
            [eigenVec,D]=eigs(K,obj.nbasis/2,'smallestabs');
%             [eigenVec2,D]=eigs(K,obj.nbasis/2);
%             eigenVec=[eigenVec eigenVec2];
            psi = K*eigenVec;

            for i = 1:size(eigenVec,2)
                b = eigenVec(:,i);
                b1 = obj.boundaryConditions.reducedToFullVector(b);
                basis{i} = reshape(b1,2,[])';
                basisVec{i}= b1;
                a = psi(:,i);
                a1=obj.boundaryConditions.reducedToFullVector(a);
                psiD{i}=reshape(a1,2,[])';
                %  bC{i} = b;

                % sF.fValues = reshape(b1,2,[])';
                % sF.mesh    = mesh;
                % bF{i} =P1Function(sF);
                % bF.plot
            end

        end

        function M = createModalPreconditioner(obj)
            for i=1:obj.nbasis
                basisVecRed(:,i)=obj.boundaryConditions.fullToReducedVector(obj.basisVec{i});
            end
            M{1} = basisVecRed';
            M{2} = inv(obj.Kmodal);
            M{3} = basisVecRed;
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
%                 function z = applyPreconditioner(obj,r)
%                     w=1;
%                     DF = obj.D-w*obj.L;
%                     DR = obj.D-w*obj.Lt;
%                     z=mldivide(DF,r);
%                     r=obj.D*z;
%                     z=mldivide(DR,r);
%                     z=w*(2-w)*z;
%                 end
        %
%                 function z = applyPreconditioner(obj,r)
%                     z=obj.D\r;
%                 end
% 
%         function z = applyPreconditioner(obj,r)
%             Lchol=obj.Lchol;
%             z = Lchol\r;
%             z = (Lchol')\z;
%         end
    end

    methods (Access = public, Static)
        function mesh = createMesh()
            % file = 'CantileverBeam_Triangle_Linear';
            % % file = 'Cantileverbeam_Quadrilateral_Bilinear';
            % a.fileName = file;
            % s = FemDataContainer(a);
            % mesh = s.mesh;


            % Generate coordinates
            x1 = linspace(0,2,20);
            x2 = linspace(0,1,20);
            % Create the grid
            [xv,yv] = meshgrid(x1,x2);
            % Triangulate the mesh to obtain coordinates and connectivities
            [F,coord] = mesh2tri(xv,yv,zeros(size(xv)),'x');

            s.coord = coord(:,1:2);
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

        function m = computeMassMatrix(mesh,displacementFun)
            s.type     = 'MassMatrix';
            s.mesh     = mesh;
            s.test      = displacementFun;
            s.trial      = displacementFun;
            s.quadratureOrder = 'QUADRATIC';
            lhs = LHSintegrator.create(s);
            m   = lhs.compute();
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

        function LHSglobal = createLHSmassGlobal(mesh,quad,material,modalFun)
            sL.material = material;
            sL.test = modalFun;
            sL.trial = modalFun;
            sL.mesh = mesh;
            sL.quadratureOrder = quad.order;
            LHSglobal = LHS_integratorMassGlobal(sL);
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

                p = r + (rsnew / rsold) * p;
                rsold = rsnew;
                iter = iter + 1;
                residual(iter) = norm(LHS*x - RHS); %Ax - b
                 err(iter)=norm(x-xsol);
                errAnorm(iter)=((x-xsol)')*LHS*(x-xsol);
                f(iter)= 0.5*(x'*LHS*x)-RHS'*x;
            end
        end



        %         function v = applyPreconditioner(M,x)
        %             for i = 1:size(M,2)
        %                 x=M{i}*x;
        %             end
        %             v = x;
        %         end


    end

end