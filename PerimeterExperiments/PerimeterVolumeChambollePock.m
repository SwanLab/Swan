classdef PerimeterVolumeChambollePock < handle

    properties (Access = private)
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
            obj.createFilterChambollePock();            
            obj.createInitialTopoloy();
            obj.obtianInitialVolum();
            obj.solveProblem();
        end

    end

    methods (Access = private)

        function init(obj)

        end

        function createMesh(obj)
            nx=60;
            ny=60;
            Lx = 1;
            Ly = 1;
            obj.mesh = TriangleMesh(Lx,Ly,nx,ny);
        end

        function createFilterChambollePock(obj)
            s.mesh  = obj.mesh;
            obj.filter = FilterPDEChambollePockContinue(s);
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

            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            s.mesh  = obj.mesh;
            filter = FilterLump(s);

            obj.chi0 = filter.compute(cFun,2);            
        end

        function obtianInitialVolum(obj)
          obj.vol0 = Integrator.compute(obj.chi0,obj.mesh,2);
        end

        function solveProblem(obj)
            z0 = obj.createSigmaFunction();
            chi = obj.chi0;
            rho = obj.chi0;
            z = z0;
            vol = obj.vol0;
            for iopt=1:4000
                chiOld = chi;
                rho = obj.proximalOfPerimeter(chi);
               % plot(rho)
                chi = obj.projectToVolumeConstraint(rho,vol);
                chi = project(chi,'P1');
                close all
               % chi.plot()
               % drawnow
                volCon = obj.computeVolumeConstraint(chi,vol);
                stopCrit = Norm(chi-chiOld,'L2')/(Norm(chiOld,'L2')+1e-5)

            end

        end
     
        function volCon = computeVolumeConstraint(obj,chi,vol)
            volCon = Integrator.compute(chi,chi.mesh,2)/vol - 1;
        end

        function rho = proximalOfPerimeter(obj,chi)
            rho = obj.filter.compute(chi);
        end

        function chi = projectToVolumeConstraint(obj,u,vol)
            chi = obj.createRho();
            lLb=0; lUb=1;
            for idic=1:100
                lam = (lLb+lUb)/2;
                chi = obj.createCharacteristicFunction(u,lam);
                constraint = obj.computeVolumeConstraint(chi,vol);
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