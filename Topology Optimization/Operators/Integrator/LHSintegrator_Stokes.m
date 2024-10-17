classdef LHSintegrator_Stokes < handle %LHSintegrator

    properties (GetAccess = public, SetAccess = private)
        M
    end

    properties (Access = private)
        dt
        mesh
        velocityFun
        pressureFun
        material
        D
    end

    methods (Access = public)

        function obj = LHSintegrator_Stokes(cParams)
            obj.init(cParams);
        end

        function LHS = compute(obj) % Venim de StokesProblem per calcular la LHS. Anem executant les funcions privades d'aquesta classe
            velLHS = obj.computeVelocityLHS(); %
            D      = obj.computeWeakDivergenceMatrix();
            prsLHS = obj.computePressureLHS(D);
            LHS = [velLHS, D; D',prsLHS];
        end

    end

%% Mètodes privats:
    methods (Access = private)
    
        function init(obj, cParams) % Li passem els paràmetres que necessitem a l'objecte d'aquesta classe (això ho fem abans de cridar al compute)
            obj.dt          = cParams.dt;
            obj.mesh        = cParams.mesh;
            obj.material    = cParams.material;
            obj.pressureFun = cParams.pressureFun;
            obj.velocityFun = cParams.velocityFun;
        end

        function LHS = computeVelocityLHS(obj) % Ens porta aquí a baix (computeVelocityLaplacian)
            K = obj.computeVelocityLaplacian();
            M = obj.computeMassMatrix();
            lhs = K + M;
            LHS = obj.symGradient(lhs);
        end

        function D = computeWeakDivergenceMatrix(obj)
            s.type = 'WeakDivergence';
            s.mesh = obj.mesh;
            s.trial = obj.pressureFun;
            s.test  = obj.velocityFun;
            LHS = LHSintegrator.create(s);
            D = LHS.compute();
        end

        function BB = computePressureLHS(obj,D)
            sz = size(D, 2);
            BB = sparse(sz,sz);
        end

        function A = symGradient(obj, B)
            A = 1/2 * (B+B');
        end

        function lhs = computeVelocityLaplacian(obj) % Laplacià
            s.type  = 'Laplacian';
            s.mesh  = obj.mesh;
            s.test  = obj.velocityFun;
            s.trial = obj.velocityFun;
            s.material = obj.material;
            LHS = LHSintegrator.create(s); % Tornem al LHSintegrator.create triant que ens doni el LHSintegrator_Laplacian. També, com així és la funció d'inicialització de
            % LHSintegrator_Laplacian, calculem també la quadratura de
            % Gauss int
            lhs = LHS.compute(); 
            lhs = obj.symGradient(lhs);
        end

        function M = computeMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = obj.velocityFun;
            s.trial = obj.velocityFun;
            s.quadratureOrder = 3;
            LHS = LHSintegrator.create(s);
            m = LHS.compute();

            dtime = obj.dt;
            M = m/dtime;
            obj.M = M;
        end

    end

end