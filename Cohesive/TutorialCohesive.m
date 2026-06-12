classdef TutorialCohesive < handle

    properties (Access = public)
        output  
        inputData
    end

    properties (Access = private)
        boundaryConditions
        functional
        solverType
        cohesiveMesh
        tractionSeparation
    end

    methods (Access = public)

        function obj = TutorialCohesive()
            obj.init();
            obj.createMesh();
            obj.defineCase();
            obj.createCohesiveFunctional()
            obj.solveCohesiveProblem()
        end

        function solveCohesiveProblem(obj)
            s.cohesiveMesh           = obj.cohesiveMesh;
            s.boundaryConditions     = obj.boundaryConditions;
            s.functional             = obj.functional;
            s.tolerance              = 1e-4;
            s.maxIter                = 100;
            s.tractionLaw            = obj.tractionSeparation;
            
            s.monitoring.set         = true;
            s.monitoring.print       = true;

            s.solverType             = obj.solverType;

            CohComp = CohesiveComputer(s);
            CohComp.compute();
            obj.output = CohComp.data;
        end

    end

    methods (Access = private)

        function init(obj)
            close all
            obj.inputData.young     = 120e20;
            obj.inputData.poisson   = 0.3;
    
            obj.inputData.lawType = 'TractionBiliniarCoupled';
            
            % TURON
            obj.inputData.tau0Normal       = 15e6;
            obj.inputData.tau0Shear        = 30e6;
            obj.inputData.firstCritEnergy  = 5*260;
            obj.inputData.secondCritEnergy = 1002;
            obj.inputData.eta              = 2;
                        
            % NO TURON
            obj.inputData.jumpCrit         = 1.25e-7;
            obj.inputData.jumpFinal        = 0.025e-3;

            % % UnitElem
            obj.inputData.Kcoh        = 1e12;
            obj.inputData.bcValues    = [0:1e-5:1.15e-4,   1.15e-4:-1e-6:0 , 0:-1e-6:-5e-5];
            obj.inputData.problemType = 'DisplacementTractionY';
            obj.inputData.l = 1;
            obj.inputData.h = 1;
            obj.inputData.yCohLine = 0;
            obj.inputData.xCohLineMax = 1000;
            obj.inputData.nx = 1;
            obj.inputData.ny = 1;

            % % DCB
            % obj.inputData.Kcoh      = 1e13;
            % obj.inputData.bcValues = [0:0.0001:1.2e-3,1.2e-3:0.00001:1.5e-3, 1.5e-3:0.00005:7e-3, 7e-3:0.00001:10e-3]; % Amb Gc = Gc
            % % obj.inputData.bcValues = [0:0.0001:2.2e-3,2.2e-3:0.00002:3.5e-3, 3.5e-3:0.00002:10e-3]; % Amb Gc = 10Gc
            % obj.inputData.problemType = 'DoubleCantileverBeam';
            % obj.inputData.l = 150e-3;
            % obj.inputData.h = 1.55*2e-3;
            % obj.inputData.yCohLine = 0.5*obj.inputData.h;
            % obj.inputData.xCohLineMax = 115e-3;
            % obj.inputData.nx = 2000;
            % obj.inputData.ny = 10;


            % % ENF
            obj.inputData.Kcoh      = 1e13;
            obj.inputData.bcValues = [ ...
                    0     : 1e-4     :6e-3,...
                    6e-3  : 2e-5     :6.5e-3,...
                    6.5e-3  : 1e-6     :9e-3];

            obj.inputData.problemType = 'EndNotchedFlex';
            obj.inputData.l = 150e-3;
            obj.inputData.h = 1.55*2e-3;
            obj.inputData.yCohLine = 0.5*obj.inputData.h;
            obj.inputData.xCohLineMax = 115e-3;
            obj.inputData.nx = 1000;
            obj.inputData.ny = 10;


    

        end
    
        function createMesh(obj)
            s.xCohLineMax = obj.inputData.xCohLineMax; s.yCohLine = obj.inputData.yCohLine;
            s.baseMesh = QuadMesh(obj.inputData.l ,obj.inputData.h, obj.inputData.nx, obj.inputData.ny);
            obj.cohesiveMesh = NewCohesiveMesh(s);
        end

        function defineCase(obj)
            obj.solverType = 'Newton';
            obj.createBoundaryConditions();
            obj.createTractionSeparation();
        end

        function createBoundaryConditions(obj)
            bc.type = obj.inputData.problemType;
            bc.values = obj.inputData.bcValues;
            obj.boundaryConditions  = BoundaryConditionsCreator(obj.cohesiveMesh.fullMesh,bc);
        end

        function createTractionSeparation(obj)
            s = obj.inputData;
            obj.tractionSeparation = CohesiveTractionSeparation(s);
        end

        function createCohesiveFunctional(obj)
            s.tractionSeparation = obj.tractionSeparation;
            s.material           = obj.createMaterial();
            s.mesh               = obj.cohesiveMesh.fullMesh;
            s.cohesiveMesh       = obj.cohesiveMesh;
            s.quadOrder          = 2;
            s.test = LagrangianFunction.create(obj.cohesiveMesh.fullMesh,2,'P1');
            obj.functional = CohesiveFunctional(s);
        end

        function m = createMaterial(obj)
            k.type    = 'ISOTROPIC';
            k.ptype   = 'ELASTIC';
            k.ndim    = obj.cohesiveMesh.fullMesh.ndim;
            k.young   = ConstantFunction.create(obj.inputData.young,obj.cohesiveMesh.fullMesh);
            k.poisson = ConstantFunction.create(obj.inputData.poisson, obj.cohesiveMesh.fullMesh);
            k.cohesiveMesh = obj.cohesiveMesh;
            m  = Material.create(k);
            m.setCohesiveMesh(obj.cohesiveMesh);
        end

        function compareDCBAnalytical(obj)

            G13 = obj.inputData.young/(2*(1 + obj.inputData.poisson));
            Gamma = 1.18*sqrt(obj.inputData.young*obj.inputData.young)/G13;
            chi = sqrt((obj.inputData.young/(11*G13))*(3 - 2*(Gamma/(1 + Gamma))^2));
            a0 = obj.inputData.l - obj.inputData.xCohLineMax;
            
            b = 1;
            
            h = 0.5*obj.inputData.h;
            Fdcb = sqrt((obj.inputData.young*b^2*h^3*obj.inputData.firstCritEnergy)/(12*(a0 + chi*h)^2));
            udcb = (8*(a0 + chi*obj.inputData.h)^3/(obj.inputData.young*b*h^3))*Fdcb;
            Cdcb = 8*(a0 + chi*h)^3/(obj.inputData.young*b*h^3);
            FAnalytical = min([0:0.01e-3:4e-3]/Cdcb,Fdcb);
            
            figure
            plot([0:0.01e-3:4e-3]*1e3,FAnalytical,'LineWidth',2)
            hold on
            grid on
            xlabel('Opening displacement [mm]')
            ylabel('Load [N]')
            title('DCB LEFM analytical solution')
            xlim([0,max(obj.inputData.bcValues)*1e3])
            ylim([0,1.05*max(FAnalytical)])
            hold on
            E  = obj.inputData.young;
            Gc = obj.inputData.firstCritEnergy;
            a0   = 35e-3;   % initial crack length [m]
            aMax = 120e-3;  % final crack length [m]
            a = linspace(a0,aMax,1000);
            % LEFM load (force per unit width)
            F = b .* sqrt(Gc .* E .* h.^3 ./ (12 .* a.^2));
            % Compliance
            C = 8 .* a.^3 ./ (E .* b .* h.^3);
            % Opening displacement
            delta = C .* F;
            plot(delta*1e3,F,'LineWidth',2)
            grid on

            Fsim = obj.output(:,1);
            Usim = obj.output(:,2)*1e3;

            plot(2*Usim,Fsim);
        end

    end

end