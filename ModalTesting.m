classdef ModalTesting < handle

    properties (Access = public)

    end

    properties (Access = private)
        nDimf
        mesh
        bMesh
        nbound
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
        RBdomain
        RBboundary
        interfaceModes
        interfaceFunTyp
        phiB
        psi
    end

    properties (Access = private)

    end

    methods (Access = public)

        function obj = ModalTesting()
            close all;
            obj.init()
            obj.mesh  = obj.createMesh();
            obj.bMesh = obj.mesh.createBoundaryMesh();
            obj.nbound = size(obj.bMesh);
            rawBc     = obj.createRawBoundaryConditions();
            obj.boundaryConditions = obj.createBoundaryConditions(rawBc);
            obj.quad = Quadrature.set(obj.mesh.type);
            obj.quad.computeQuadrature('QUADRATIC');
            obj.material = obj.createMaterial();

            dispFun = P1Function.create(obj.mesh, obj.nDimf);

            K    = obj.computeStiffnessMatrix(obj.mesh,obj.material,dispFun);
            obj.Kred = obj.boundaryConditions.fullToReducedMatrix(K);
            %             R = sprand((obj.Kred));
            %            obj.Kred=0.5*(R'*R);
            obj.D = diag(diag(obj.Kred));
            obj.L = tril(obj.Kred,-1);
            obj.Lt= obj.L';
%             obj.Lchol=ichol(obj.Kred);

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

%             [xCG,residualCG,errCG,errACG]   = obj.conjugateGradient(obj.Kred,Fred,x);

%             CG.x        = xCG;
%             CG.residual = residualCG;
%             CG.err      = errCG;
%             CG.errA     = errACG;

%             [xPCG,residualPCG,errPCG,errAPCG] = obj.preconditionedConjugateGradient(obj.Kred,Fred,x);

%             PCGbasis20.x        = xPCG;
%             PCGbasis20.residual = residualPCG;
%             PCGbasis20.err      = errPCG;
%             PCGbasis20.errA     = errAPCG;

            obj.computeBasisRom();
            [Hhat,H,G,T,That] = computeEifemMat(obj);
            nphi = obj.phiB{1}.nbasis;
            npsi = obj.psi{1}.nbasis;
            nrb  = obj.RBboundary{1}.nbasis;
            A = [ obj.Kmodal        zeros(nphi,nrb)     -Hhat              -H               ;
                  zeros(nrb,nphi)   zeros(nrb,nrb)      -G                 zeros(nrb,npsi)  ;
                  Hhat'             G'                  zeros(nrb,nrb)     zeros(nrb,npsi)  ;
                  H'                zeros(npsi,nrb)     zeros(npsi,nrb)    zeros(npsi,npsi)];
                
            b    = zeros(size(A(:,1),1),obj.interfaceModes{1}.nbasis);
            rv   = That*eye(obj.interfaceModes{1}.nbasis,obj.interfaceModes{1}.nbasis);
            psiv = T*eye(obj.interfaceModes{1}.nbasis,obj.interfaceModes{1}.nbasis);
            b(nphi+nrb+1:end,:) = [rv;psiv];
            x = A\b;
            kcoarse = x(nphi+nrb+1:end,:);
            basis   = cell2mat(obj.basisVec);
            Udef    = basis*inv(H)*T;
            for i = 1:3
                p1test=obj.RBdomain.basisFunctions{i}.project('P1');
                RB(:,i)=reshape(p1test.fValues',[],1);
            end
%              Urb     = obj.RBdomain*inv(G)*(That-Hhat*inv(H)*T);
            Urb     = RB*inv(G)*(That-Hhat'*inv(H)*T);
            kcoarse = (Udef')*K*Udef;

            rom.kcoarse = kcoarse;
            rom.Udef    = Udef;
            rom.Urb     = Urb;
            rom.ndimf   = dispFun.ndimf;
            rom.mesh    = obj.mesh;
%             rom.basis   = basis;

            save("testEIFEM1.mat","rom")

            %             ModalTesting.plotRes(residualCG,residualPCG,errCG,errPCG,errACG,errAPCG)



        end

    end

    methods (Access = private)

        function init(obj)
            obj.nDimf   = 2;
            obj.nbasis = 5;
            funTyp='P1';
            for i=1:obj.nbasis
                obj.functionType{i}=funTyp;
            end
            intFunTyp='P1';
            for i=1:4*obj.nDimf
                obj.interfaceFunTyp{i}=intFunTyp;
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
            bc.pointload(:,3) = -1/length(nodes(forceNodes));
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
            [eigenVec,D]=eigs(K,obj.nbasis,'smallestabs');
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



        function computeBasisRom(obj)
            obj.RBdomain       = obj.createRBfun(obj.mesh);
            obj.RBboundary     = obj.createRBfun(obj.bMesh);
            obj.phiB           = obj.createBoundaryDefFun();
            obj.psi            = obj.phiB;
            obj.interfaceModes = obj.computeBasisInterface();   
        end

        function modes = computeBasisInterface(obj)
            cCoord = obj.getFictitiousCornerCoord();
            lx = max(cCoord(:,1))-min(cCoord(:,1));
            ly = max(cCoord(:,2))-min(cCoord(:,2));
            cCoord = [ cCoord(end,:); cCoord];
            for ibound=1:size(obj.bMesh,1)
                ibasis=1;
                for icorner=2:size(cCoord,1)
                    x = obj.bMesh{ibound}.mesh.coord;
                    for idim=1:obj.nDimf
                        f=zeros(size(x));

                        if x(2,1) == cCoord(icorner,1)
                            f(:,idim) = 1-(abs((x(:,2)-cCoord(icorner,2)))/ly);
                            basis{ibasis}=f;
                            ibasis=ibasis+1;
                            %                            imodes{ibound}{icorner,idim}= ModalFunction.create(obj.bMesh{ibound}.mesh,f,obj.interfaceFunTyp);
                            %                             f(:,2) = zeros(size(x(:,1)));
                        elseif x(2,2) == cCoord(icorner,2)
                            %                             f(:,2) = zeros(size(x(:,1)));
                            f(:,idim) = 1-(abs((x(:,1)-cCoord(icorner,1)))/lx);
                            basis{ibasis}=f;
                            ibasis=ibasis+1;
                            %                             imodes{ibound}{icorner,idim}= ModalFunction.create(obj.bMesh{ibound}.mesh,f,obj.interfaceFunTyp);
                        else
                            basis{ibasis}=f;
                            ibasis=ibasis+1;
                            %                             imodes{ibound}{icorner,idim}= ModalFunction.create(obj.bMesh{ibound}.mesh,f,obj.interfaceFunTyp);

                        end
                    end

                end
                modes{ibound}= ModalFunction.create(obj.bMesh{ibound}.mesh,basis,obj.interfaceFunTyp);
            end

        end

        function coord = getFictitiousCornerCoord(obj)
            xmax = max(obj.mesh.coord(:,1));
            xmin = min(obj.mesh.coord(:,1));
            ymax = max(obj.mesh.coord(:,2));
            ymin = min(obj.mesh.coord(:,2));
            coord(1,1) = xmax;
            coord(1,2) = ymin;
            coord(2,1) = xmax;
            coord(2,2) = ymax;
            coord(3,1) = xmin;
            coord(3,2) = ymax;
            coord(4,1) = xmin;
            coord(4,2) = ymin;

        end

        function defFun = createBoundaryDefFun(obj)
            for imesh=1:obj.nbound
                nodes = unique(obj.bMesh{imesh}.globalConnec);

                for ibasis=1:length(obj.basisFvalues)
                    basisb{ibasis} = obj.basisFvalues{ibasis}(nodes,:);
                end

                meshb = obj.bMesh{imesh}.mesh;
                defFun{imesh} = ModalFunction.create(meshb,basisb,obj.functionType);
            end
        end

        function [Hhat,H,G,T,That] = computeEifemMat(obj)
                Hhat = obj.computeHhat();
                H    = obj.computeH();
                G    = obj.computeG();
                T    = obj.computeT();
                That = obj.computeThat();
                 
        end

        function Hhat = computeHhat(obj)
            Hhat = zeros(obj.phiB{1}.nbasis,obj.RBboundary{1}.nbasis);
            for ibound=1:obj.nbound
                sL.test  = obj.phiB{ibound};
                sL.trial = obj.RBboundary{ibound};
                sL.mesh  = obj.bMesh{ibound}.mesh;
                sL.quadratureOrder = obj.quad.order;
                lhs2{ibound} = LHS_integratorMassGlobal(sL);
                Hhat = Hhat + lhs2{ibound}.compute();
%                 Hhat{ibound} = lhs2{ibound}.compute();
            end
        end

        function H = computeH(obj)
            H = zeros(obj.psi{1}.nbasis,obj.phiB{1}.nbasis);
            for ibound=1:obj.nbound
                sL.test  = obj.psi{ibound};
                sL.trial = obj.phiB{ibound};
                sL.mesh  = obj.bMesh{ibound}.mesh;
                sL.quadratureOrder = obj.quad.order;
                lhs2{ibound} = LHS_integratorMassGlobal(sL);
                H = H + lhs2{ibound}.compute();
%                 H{ibound} = lhs2{ibound}.compute();
            end
        end

        function G = computeG(obj)
            G = zeros(obj.RBboundary{1}.nbasis,obj.RBboundary{1}.nbasis);
            for ibound=1:obj.nbound
                sL.test  = obj.RBboundary{ibound};
                sL.trial = obj.RBboundary{ibound};
                sL.mesh  = obj.bMesh{ibound}.mesh;
                sL.quadratureOrder = obj.quad.order;
                lhs2{ibound} = LHS_integratorMassGlobal(sL);
                G = G + lhs2{ibound}.compute();
%                 G{ibound} = lhs2{ibound}.compute();
            end
        end

        function T = computeT(obj)
            T = zeros(obj.psi{1}.nbasis,obj.interfaceModes{1}.nbasis);
            for ibound=1:obj.nbound
                sL.test  = obj.psi{ibound};
                sL.trial = obj.interfaceModes{ibound};
                sL.mesh  = obj.bMesh{ibound}.mesh;
                sL.quadratureOrder = obj.quad.order;
                lhs2{ibound} = LHS_integratorMassGlobal(sL);
                T = T + lhs2{ibound}.compute();
%                 T{ibound} = lhs2{ibound}.compute();
            end
        end

        function That = computeThat(obj)
            That = zeros(obj.RBboundary{1}.nbasis,obj.interfaceModes{1}.nbasis);
            for ibound=1:obj.nbound
                sL.test  = obj.RBboundary{ibound};
                sL.trial = obj.interfaceModes{ibound};
                sL.mesh  = obj.bMesh{ibound}.mesh;
                sL.quadratureOrder = obj.quad.order;
                lhs2{ibound} = LHS_integratorMassGlobal(sL);
                That = That + lhs2{ibound}.compute();
%                 T{ibound} = lhs2{ibound}.compute();
            end
        end


    end

    methods (Access = public, Static)
        function mesh = createMesh()
            % file = 'CantileverBeam_Triangle_Linear';
            % % file = 'Cantileverbeam_Quadrilateral_Bilinear';
            % a.fileName = file;
            % s = FemDataContainer(a);
            % mesh = s.mesh;

%             filename   = 'lattice_ex1';
%             a.fileName = filename;
%             femD       = FemDataContainer(a);
%             mesh       = femD.mesh;
%             bS         = mS.createBoundaryMesh();


            % Generate coordinates
            x1 = linspace(0,2,50);
            x2 = linspace(0,1,50);
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

        function plotRes(residualCG,residualPCG,errCG,errPCG,errACG,errAPCG)
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

        function rbFun = createRBfun(mesh)
            %domain
            if size(mesh,1)==1
                centroid = [sum(mesh.coord(:,1))/mesh.nnodes,sum(mesh.coord(:,2))/mesh.nnodes];
                rbFun = RigidBodyFunction.create(mesh,centroid);
            else
                %boundary
                for i=1:size(mesh,1)
                    meshb = mesh{i}.mesh;
                    centroid = [sum(meshb.coord(:,1))/meshb.nnodes,sum(meshb.coord(:,2))/meshb.nnodes];
                    rbFun{i} = RigidBodyFunction.create(meshb,centroid);
                end
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