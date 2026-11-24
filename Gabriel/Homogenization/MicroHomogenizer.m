classdef MicroHomogenizer < handle

    properties (Access = private)
        E
        nu
        voidRatio
    end

    properties (Constant, Access = private)
        pointTolerance = 1e-12
    end

    methods (Access = public)

        function obj = MicroHomogenizer(cParams)
            obj.E         = cParams.E;
            obj.nu        = cParams.nu;
            obj.voidRatio = cParams.voidRatio;
        end

        function C = compute(obj,periodicCell,dens)
            mesh     = periodicCell.mesh;
            material = obj.createMaterial(mesh);
            s.mesh               = mesh;
            s.material           = [];
            s.scale              = 'MICRO';
            s.dim                = '2D';
            s.boundaryConditions = obj.createBoundaryConditions(periodicCell);
            s.solverCase         = DirectSolver();
            s.solverType         = 'REDUCED';
            s.solverMode         = 'FLUC';
            fem = ElasticProblemMicro(s);
            fem.updateMaterial(material(dens));
            fem.solve();
            C = fem.Chomog/mesh.computeVolume();
        end

    end

    methods (Access = private)

        function C = createMaterial(obj,mesh)
            N            = mesh.ndim;
            [muA,kappaA] = obj.computeModuli(obj.voidRatio*obj.E,N);
            [muB,kappaB] = obj.computeModuli(obj.E,N);
            mu     = @(rho) Expand(SimpAllInterpolator.computeMu(muA,muB,kappaA,kappaB,rho,N),4);
            kappa  = @(rho) Expand(SimpAllInterpolator.computeKappa(muA,muB,kappaA,kappaB,rho,N),4);
            lambda = @(rho) kappa(rho) - (2/N)*mu(rho);
            I      = ConstantFunction.create(eye4D(N),mesh);
            IxI    = ConstantFunction.create(kronEye(N),mesh);
            C      = @(rho) 2*mu(rho).*I + lambda(rho).*IxI;
        end

        function [mu,kappa] = computeModuli(obj,E,N)
            mu    = E/(2*(1+obj.nu));
            kappa = E/(N*(1-(N-1)*obj.nu));
        end

        function bc = createBoundaryConditions(obj,periodicCell)
            mesh             = periodicCell.mesh;
            vertexCoord      = mesh.coord(periodicCell.vertexNodes,:);
            sDir.domain      = @(coor) obj.isAtPoints(coor,vertexCoord);
            sDir.direction   = 1:mesh.ndim;
            sDir.value       = 0;
            s.dirichletFun   = DirichletCondition(mesh,sDir);
            s.pointloadFun   = [];
            s.periodicFun    = 1;
            s.masterSlave    = periodicCell.masterSlave;
            s.mesh           = mesh;
            bc = BoundaryConditions(s);
        end

        function is = isAtPoints(obj,coor,points)
            is = false(size(coor,1),1);
            for i = 1:size(points,1)
                is = is | vecnorm(coor-points(i,:),2,2) < obj.pointTolerance;
            end
        end

    end

end