classdef PerimeterVolumeChambollePock < handle

    properties (Access = private)
        eps
        thetaRel
        tauG
        tauF
        mesh
        filter
        chi0
        vol0
    end

    properties (Access = private)

    end

    methods (Access = public)

        function obj = PerimeterVolumeChambollePock()
            obj.init()
            k=[1,1]; k=k'/norm(k); kx=k(1); ky=k(2);
            kperp=[-ky;kx];            
            obj.createMesh();
            obj.createFilter();
            obj.createInitialTopoloy();
            obj.obtianInitialVolum();
            obj.createProximals();
            obj.solveProblem();
        end

    end

    methods (Access = private)

        function init(obj)
            obj.tauF = 0.0025;
            obj.tauG = 0.25;
            obj.thetaRel = 1;
            obj.eps  = 1;
        end

        function createMesh(obj)
            nx=50;
            ny=50;
            Lx = 1;
            Ly = 1;
            obj.mesh = TriangleMesh(Lx,Ly,nx,ny);
        end

        function createFilter(obj)
            s.trial = obj.createRho();
            s.mesh  = obj.mesh;
            obj.filter = FilterLump(s);
        end

        function g = createGeometricalFunction(obj)
            fracx = 1/3;
            fracy = 1/3;
            Lx = max(obj.mesh.coord(:,1));
            Ly = max(obj.mesh.coord(:,2));
            cx = Lx/2;
            cy = Ly/2;
            w  = fracx*Lx;
            h = fracy*Ly;
            s.xSide = w;
            s.ySide = h;
            s.xCoorCenter = cx;
            s.yCoorCenter = cy;
            s.type = 'Rectangle';
            g = GeometricalFunction(s);
        end

        function createInitialTopoloy(obj)
            g = obj.createGeometricalFunction;
            ls = g.computeLevelSetFunction(obj.mesh);
            sU.backgroundMesh = obj.mesh;
            sU.boundaryMesh   = obj.mesh.createBoundaryMesh;
            uMesh             = UnfittedMesh(sU);
            uMesh.compute(ls.fValues);
            cFun = CharacteristicFunction.create(uMesh);
            obj.chi0 = obj.filter.compute(cFun,2);            
        end

        function obtianInitialVolum(obj)
          obj.vol0 = Integrator.compute(obj.chi0,obj.mesh,2);
        end

        function [proxF,proxG] = createProximals(obj)
            %proxF = @(z)  proximalDroplet(z,tauF,k,alpha,ep);
            epsF    = 1; 
            proxF   = @(z)  obj.proximalEllipse(z,obj.tauF,epsF,obj.mesh);           
            tauGEps = obj.tauG/obj.eps^2;
            proxG   = @(rho,chi) obj.proximalL2Projection(rho,chi,tauGEps);
        end

        function solveProblem(obj)
            [proxF,proxG] = obj.createProximals();
            z0 = obj.createSigmaFunction();
            chi = obj.chi0;
            rho = obj.chi0;
            z = z0;
            for iopt=1:20
                [rho,z] = obj.proximalOfPerimeter(chi,rho,z0,proxF,proxG,obj.tauF,obj.tauG,obj.thetaRel);
                plot(rho)
                chi = obj.projectToVolumeConstraint(rho,vol);
                close all
                chi.plot()
                drawnow
                volCon = obj.computeVolumeConstraint(chi,vol);
            end

        end
     
        function volCon = computeVolumeConstraint(chi,vol)
            volCon = Integrator.compute(chi,chi.mesh,2)/vol - 1;
        end

        function s = proximalEllipse(obj,z,tau,alpha,m)
            A = [1 0; 0 1];
            I = eye(2);
            r = alpha^2/tau;
            dm = obj.createTensorFunction(A+r*I);
            s = obj.createSigmaFunction();
            M  = IntegrateLHS(@(u,v) DP(v,DP(dm,u)),s,s,m,'Domain');
            %LHS = diag(sum(M));
            LHS = M;
            F  = IntegrateRHS(@(v) DP(v,z),s,m,'Domain');
            sV = full(LHS\F);
            sV = reshape(sV,[s.ndimf,m.nnodes])';
            s.setFValues(sV);

        end

        function s = proximalDroplet(obj,z,tau,k,alpha,ep)
            s1 = z/(1+tau);
            zk  = z*k(:);
            z2 = sum(z.^2,2);
            tA   = tau/(alpha^2*(1+tau));
            tB   = sqrt((alpha^2-1)./(z2-zk.^2+ep));
            deltaS = tA*(z + (alpha^2-2)*zk.*k.' - tB.*((z2-2*zk.^2).*k.' + zk.*z));
            s      = s1 + (alpha*zk - sqrt(z2) > 0).*deltaS;
        end

        function rhoN = proximalL2Projection(obj,rho,Chi,tauG)
            rhoN = obj.filter.compute((1/(1+tauG))*(rho + tauG*Chi),2);
        end

        function [u,z] = solveWithChambollePockAlgorithm(obj,u0,z0,proxF,proxG,tauF,tauG,thetaRel)
            u  = u0;
            uN = u0;
            z   = z0;
            for kcp=1:100
                z      = proxF(z + tauF.*Grad(uN));
                uOld   = u;
                u      = proxG(u - tauG*Divergence(z));
                uN     = project(u + thetaRel.*(u - uOld),'P1');
            end
        end

        function [u,z] = proximalOfPerimeter(obj,chi,u0,z0,proxF,proxGX,tauF,tauG,thetaRel)
            proxG = @(rho) proxGX(rho,chi);
            [u,z] = obj.solveWithChambollePockAlgorithm(u0,z0,proxF,proxG,tauF,tauG,thetaRel);
        end

        function Af = createTensorFunction(obj,A)
            op = @(xV) obj.evaluate(xV,A);
            Af = DomainFunction.create(op,obj.mesh,2);
        end

        function fV = evaluate(obj,xV,A)
            nGauss = size(xV,2);
            nElem  = obj.mesh.nelem;
            fV = repmat(A,[1 1 nGauss nElem]);
        end


        function chi = projectToVolumeConstraint(obj,u,vol)
            chi = createRho(u.mesh);
            lLb=0; lUb=1;
            for idic=1:100
                lam = (lLb+lUb)/2;
                chi = createCharacteristicFunction(u,lam);
                constraint = computeVolumeConstraint(chi,vol);
                if (constraint > 0)
                    lLb =lam;
                else
                    lUb =lam;
                end
            end
        end


        function chi = createCharacteristicFunction(obj,u,lambda)
            op = @(xV) obj.evaluateCharact(xV,u,lambda);
            chi = DomainFunction.create(op,u.mesh,1);
        end

        function fV = evaluateCharact(obj,xV,u,lambda)
            fV = real(u.evaluate(xV)-lambda>0);
        end

       function s = createSigmaFunction(obj)
            s = LagrangianFunction.create(obj.mesh,2,'P1');
        end

        function rho = createRho(obj)
            rho = LagrangianFunction.create(obj.mesh,1,'P1');
        end


        end

end