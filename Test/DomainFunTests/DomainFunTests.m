classdef DomainFunTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        cases = {'DDP','Grad','Partial','SymGrad','Divergence','VolumetricElas','DeviatoricElas'}
    end

    methods (Test, TestTags = {'DomainFun'})
        function test2D(testCase, cases)
            filename = ['testDomainFun',cases,'2D'];
            m        = testCase.obtain2DTestMesh();
            E        = 1;
            nu       = 1/3;
            C        = testCase.computeMaterial(m,E,nu);
            bc       = testCase.createBC(m,'2D');
            fem      = testCase.solveElasticProblem(m,C,'2D',bc);
            kappa    = C.computeKappaFromYoungAndPoisson(E,nu,2);
            kappaFun = ConstantFunction.create(kappa,m);
            mu       = C.computeMuFromYoungAndPoisson(E,nu);
            muFun    = ConstantFunction.create(mu,m);
            switch cases
                case 'DDP'
                    domainFun = DDP(C,SymGrad(fem.uFun));
                case 'Grad'
                    domainFun = Grad(fem.uFun);
                case 'Partial'
                    domainFun = Partial(fem.uFun,1);
                case 'SymGrad'
                    domainFun = SymGrad(fem.uFun);
                case 'Divergence'
                    domainFun = Divergence(fem.uFun);
                case 'VolumetricElas'
                    domainFun = VolumetricElasticEnergyDensity(fem.uFun,kappaFun);
                case 'DeviatoricElas'
                    domainFun = DeviatoricElasticEnergyDensity(fem.uFun,muFun);
            end
            quad = Quadrature.create(m, 5);
            xV   = quad.posgp;
            xNew = domainFun.evaluate(xV);
            load(filename,'xRef');
            err  = norm(xNew(:)-xRef(:))/norm(xRef(:));
            tol  = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

        function test3D(testCase, cases)
            filename = ['testDomainFun',cases,'3D'];
            m        = testCase.obtain3DTestMesh();
            E        = 1;
            nu       = 1/3;
            C        = testCase.computeMaterial(m,E,nu);
            bc       = testCase.createBC(m,'3D');
            fem      = testCase.solveElasticProblem(m,C,'3D',bc);
            kappa    = C.computeKappaFromYoungAndPoisson(E,nu,3);
            kappaFun = ConstantFunction.create(kappa,m);
            mu       = C.computeMuFromYoungAndPoisson(E,nu);
            muFun    = ConstantFunction.create(mu,m);
            switch cases
                case 'DDP'
                    domainFun = DDP(C,SymGrad(fem.uFun));
                case 'Grad'
                    domainFun = Grad(fem.uFun);
                case 'Partial'
                    domainFun = Partial(fem.uFun,1);
                case 'SymGrad'
                    domainFun = SymGrad(fem.uFun);
                case 'Divergence'
                    domainFun = Divergence(fem.uFun);
                case 'VolumetricElas'
                    domainFun = VolumetricElasticEnergyDensity(fem.uFun,kappaFun);
                case 'DeviatoricElas'
                    domainFun = DeviatoricElasticEnergyDensity(fem.uFun,muFun);
            end
            quad = Quadrature.create(m, 5);
            xV   = quad.posgp;
            xNew = domainFun.evaluate(xV);
            load(filename,'xRef');
            err  = norm(xNew(:)-xRef(:))/norm(xRef(:));
            tol  = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end
    end

    methods (Static, Access = private)
        function m = obtain2DTestMesh()
            m = QuadMesh(2,1,20,10);
        end

        function m = obtain3DTestMesh()
            m = HexaMesh(2,1,1,10,5,5);
        end

        function mat = computeMaterial(m,E1,nu1)
            E   = ConstantFunction.create(E1,m);
            nu  = ConstantFunction.create(nu1,m);
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = m.ndim;
            s.young   = E;
            s.poisson = nu;
            mat       = Material.create(s);
        end

        function bc = createBC(m,dim)
            xMax    = max(m.coord(:,1));
            isDir   = @(coor)  abs(coor(:,1))==0;
            isForce = @(coor)  abs(coor(:,1))==xMax;

            sDir{1}.domain    = @(coor) isDir(coor);
            switch dim
                case '2D'
                    sDir{1}.direction = [1,2];
                case '3D'
                    sDir{1}.direction = [1,2,3];
            end
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isForce(coor);
            sPL{1}.direction = 2;
            sPL{1}.value     = -1;

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(m, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;

            pointloadFun = [];
            for i = 1:numel(sPL)
                pl = PointLoad(m, sPL{i});
                pointloadFun = [pointloadFun, pl];
            end
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh = m;
            bc = BoundaryConditions(s);
        end

        function fem = solveElasticProblem(m,mat,dim,bc)
            s.mesh = m;
            s.scale = 'MACRO';
            s.material = mat;
            s.dim = dim;
            s.boundaryConditions = bc;
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'DIRECT';
            fem = ElasticProblem(s);
            fem.solve();
        end

    end

end