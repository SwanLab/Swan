classdef mainTFGMultimaterial < handle

    properties (Access = public)
        psi
        displacements
        forces
        designVariable
        tensor
        strain
    end

    properties (Access = private)
        dC
        compliance
        dJCompliance
        tgamma
        mesh
        meshSwan
        mat
        params
        bc
        bcSwan
        pdeCoeff
        initialLevelSet
        pfree
        volume
        charFunc
        optParams
        unitM
        shapeFunc
        Epot
        DT
        nMat
        cosAngle
        thetaAngle
        k
        fiFunction
        phold
        iter
        area
        unitMSeba
        displFun
        monitoring
    end


    methods (Access = public)
    
        function obj = mainTFGMultimaterial()
            obj.init();
            obj.uploadParameters();
            obj.createMesh();
            obj.createBoundaryConditions();
            obj.createMaterial();
            obj.computeInitialLevelSet();
            obj.plotMesh();
            obj.computeFreeNodes();
            obj.solveInitialElasticProblem();
            obj.createOptimizationParameters();
            obj.createMassMatrix();
            %obj.createMassMatrixSwan();
            obj.compute();
        end

        function C = computeElasticTensor(obj)
            s.psi               = obj.psi;
           % s.mesh              = obj.mesh;
            s.matProp           = obj.mat;
            s.pdeCoeff          = obj.pdeCoeff;
            s.bc                = obj.bcSwan;
            s.m                 = obj.meshSwan;
            s.designVariable    = obj.designVariable;
            
            constituitiveTensor = ElasticTensorComputer(s);
            C = constituitiveTensor.C;
        end

        function [U,F] = solveFEM(obj)
            s.mesh               = obj.meshSwan;
            s.scale              = 'MACRO';
            s.material           = obj.tensor;
            s.dim                = '2D';
            s.boundaryConditions = obj.bcSwan;
            s.solverType         = 'REDUCED';
            s.solverMode         = 'DISP';

            fem         = ElasticProblem(s);
            fem.solve();
            displ      = fem.uFun.fValues; 
            obj.displFun = fem.uFun;
            U          = reshape(displ, [size(displ,1)*2, 1]);
            force      = fem.forces; 
            Fx         = force(1:2:end); % Fx values are at odd indices
            Fy         = force(2:2:end); % Fy values are at even indices
            forcesVect = [Fx Fy];
            F          = reshape(forcesVect, [size(displ,1)*2, 1]);
        end

        function volume = computeVolume(obj)
            [~,tfi] = computeCharacteristicFunction(obj);  
            s.tfi   = tfi;
            s.mesh  = obj.meshSwan;
            s.area = obj.area;
            vol     = VolumeComputer(s);
            volume  = vol.computeVolume();     
        end

        function [fi, tfi] = computeCharacteristicFunction(obj)
            s.psi            = obj.psi;
            s.p              = obj.mesh.p;
            s.t              = obj.mesh.t;
            s.designVariable = obj.designVariable;
            s.m              = obj.meshSwan;
    
            charfun          = CharacteristicFunctionComputer(s); 
            [fi,tfi]         = charfun.computeFiandTfi();
        end

        function [sf, energyPot] = updateShapeFunctions(obj)
            U              = obj.displacements;
            F              = obj.forces;
            V              = obj.volume;
            parameters     = obj.params;
            TOParams       = obj.optParams;
            
            [sf,energyPot] = shfunc(F,U,V,parameters,TOParams);
        end

        function dt = computeTopologicalDerivative(obj)
            s.meshSeba   = obj.mesh;
            s.U          = obj.displacements;
            s.volume     = obj.volume;
            s.mat    = obj.mat;
            s.psi         = obj.psi;
            s.designVariable    = obj.designVariable;
            s.mesh          = obj.meshSwan;
            s.penalization = obj.params.penalization;
            s.penalty = obj.params.penalty;
            s.volfrac = obj.params.volfrac;
            s.auglag = obj.params.auglag;
            s.max_vol = obj.optParams.max_vol; 
            s.energy0 = obj.optParams.energy0;
            s.nMat = obj.nMat;
            s.mat = obj.mat;
            s.psi = obj.psi;     
            s.designVariable = obj.designVariable;

            topder = TopologicalDerivativeComputer(s);
            obj.dC = topder.dC;
            obj.strain = topder.strain;
            obj.tgamma = topder.tgamma;
            %obj.computeStrain();
            obj.computeComplianceAndGradient();
            TD = obj.computeVolumeConstraintinDT();
            dt = obj.smoothTopologicalDerivative(TD);
            
            dt = dt/normL2(obj.unitM,dt);               
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
            obj.nMat = 4;
        end

        function compute(obj)
            
            obj.nMat = 4;
            p = obj.mesh.p;

            ls = ones(size(p,2),(obj.nMat-1)); 
            ls(:,1) = -1;

            obj.psi = ls;
            
            obj.updateDesignVariable();
            %designVar = ls.designVariable;
            
            vol = computeVolume(obj);
            obj.volume = vol;

            obj.k = 1; obj.iter = 0;
            gsf = []; gEpot = []; gth = []; gvol = [];
            
            % Shape functional and potential energy
            [sf, energyPot] = updateShapeFunctions(obj);
            obj.shapeFunc = sf;
            obj.Epot = energyPot;

            % Compute topological derivative
            obj.DT = computeTopologicalDerivative(obj);
            obj.DT(obj.phold,:) = obj.psi(obj.phold,:); % freeze dt function

            % Compute cosinus and theta
            obj.computeThetaAndCosinus;
            theta = obj.thetaAngle;

            % Compute index control 
            difvol = obj.volume(1:end-1)-obj.optParams.voltarget; 
            ic = (abs(difvol) > obj.optParams.volstop); %index control
             
           % UPDATE MONITORING ITER 0
            obj.createMonitoring();

            % We start the loop    
            while obj.hasNotFinished(theta,ic)    
            % while and( or(any(ic),theta > params.stop) , k/2 > params.kmin)           
                
                obj.iter = obj.iter + 1;

                % First update psi    
                obj.updateDesignVariable();
            
                % Then compute again the elastic tensor 
                obj.tensor = computeElasticTensor(obj);

                % Solve elastic problem
                [U,F] = solveFEM(obj);
                obj.displacements = U;
                obj.forces = F;

                % Update volume, energy and shape functionals
                vol = computeVolume(obj);
                obj.volume = vol;
                [sf, energyPot] = updateShapeFunctions(obj);
                obj.shapeFunc = sf;
                obj.Epot = energyPot;

                % Update topological derivative
                obj.DT = computeTopologicalDerivative(obj);
                obj.DT(obj.phold,:) = obj.psi(obj.phold,:); % freeze dt function

                % Update cosinus and theta
                obj.computeThetaAndCosinus();
                theta = obj.thetaAngle;

                % Perform a line-search
                sfold = sf; psiold = obj.psi; sf = sf + 1; obj.k = min(1-(1e-5),1.5*obj.k); % trick
                     
                while and((sf > sfold) , obj.k>obj.params.kmin)

                    % Update level set using slerp
                    theta = obj.thetaAngle;
                    free = obj.pfree;
                    dt = obj.DT;
                    s.mesh = obj.meshSwan;
                    s.tau = obj.k;
                    
                    primalUpdater = SLERP(s);
                    obj.psi = primalUpdater.update(dt,psiold);

                    % obj.psi(free,:)  = (sin((1-obj.k)*theta)*psiold(free,:)...
                    %         +  sin(obj.k*theta)*dt(free,:))./sin(theta);
                        
                    % Solve elastic problem again    
                    obj.updateDesignVariable();
                
                    obj.tensor = computeElasticTensor(obj);
    
                    [U,F] = solveFEM(obj);
                    obj.displacements = U;
                    obj.forces = F;

                     % Update volume, energy and shape functionals
                    vol = computeVolume(obj);
                    obj.volume = vol;
                    [sf, energyPot] = updateShapeFunctions(obj);
                    obj.shapeFunc = sf;
                    obj.Epot = energyPot;

                    obj.k = obj.k / 2;
                end
                    %psi = psi/normL2(unitM,psi);
                    obj.k = obj.k * 2;

                % Update some parameters
                gsf  = [gsf,sf];
                gEpot = [gEpot,obj.Epot];
                gth  = [gth,obj.thetaAngle];
                gvol = [gvol,obj.volume(1:end-1)'*100/obj.optParams.max_vol];

                % Compute characteristic function for plot
                obj.updateDesignVariable();

                [fi, ~] = computeCharacteristicFunction(obj);
                obj.fiFunction = fi';

                

                % Plot the evolution of iterations
                obj.updatePlotAndDisplay();
                
                %Update augmented lagrangian parameters
                obj.increaseLagrangianMultiplier();
                
                % UPDATE MONITORING CURRENT ITERATION
                obj.updateMonitoring();

            end
        end

        function isNotFinished = hasNotFinished(obj,theta,ic)
            isNotFinished = and(and( or(any(ic),theta > obj.params.stop) , obj.k/2 > obj.params.kmin), obj.iter<=4); % remove iter<=4
        end

        function uploadParameters(obj)
            parameters = ParametersComputer();
            obj.params = parameters.params;
        end

        function createMesh(obj) 
            obj.mesh     = MeshComputer();
            s.connec     = obj.mesh.t';
            s.connec     = s.connec(:,1:3);
            s.coord      = obj.mesh.p';
            
            obj.meshSwan = Mesh.create(s);
            
            % Per fer altres exemples:
            % obj.meshSwan = TriangleMesh(6,1,150,25); % Bridge
            % obj.meshSwan = TriangleMesh(2,1,100,50); % Beam
            % obj.meshSwan = TriangleMesh(2,1,100,50); % Arch

            p = obj.meshSwan.coord';
            t = obj.meshSwan.connec';
            obj.area = pdetrg(p,t);

        end

        function createBoundaryConditions(obj)
            % s.m        = obj.meshSwan; 
            % s.mesh     = obj.mesh;
            % s.g        = obj.mesh.g;
            % bounCon    = BoundaryConditionsComputer(s);
            % obj.bc     = bounCon.bc;
            % s.mesh     = obj.meshSwan; 

            s.mesh = obj.meshSwan;
            BoundCond  = BoundaryConditionsSwan(s);
            obj.bcSwan = BoundCond.createBoundaryConditions();

             % Per fer altres exemples:
            %obj.bcSwan = BoundCond.createBoundaryConditionsTutorialBeam();
            %obj.bcSwan = BoundCond.createBoundaryConditionsTutorialBridge();
            %obj.bcSwan = BoundCond.createBoundaryConditionsTutorialArch();
        end

        function createMaterial(obj)
            matProp      = MaterialPropertiesComputer(); 
            obj.mat.A    = matProp.matA; % MATERIAL 1
            obj.mat.B    = matProp.matB; % MATERIAL 2
            obj.mat.C    = matProp.matC; % MATERIAL 3
            obj.mat.D    = matProp.matD; % VOID
            s.mat        = obj.mat;
            s.m          = obj.meshSwan;
            
            obj.pdeCoeff = PDECoefficientsComputer(s);  
        end

        function computeInitialLevelSet(obj)
            s.mesh                 = obj.meshSwan;
            s.type                 = 'Full';
            lsFun{1}               = -ones(size(s.mesh.coord,1),1);
            lsFun{2}               = ones(size(s.mesh.coord,1),1);
            lsFun{3}               = ones(size(s.mesh.coord,1),1);
            
            s.type                 = 'LevelSet';
            s.plotting             = false;
            s.fValues              = lsFun{1};
            s.order                = 'P1';
            s.fun                  = LagrangianFunction(s);
            ls1                    = DesignVariable.create(s);
            obj.designVariable{1} = ls1;
            
            s.fValues              = lsFun{2};
            s.fun                  = LagrangianFunction(s);
            ls1                    = DesignVariable.create(s);
            obj.designVariable{2} = ls1;
            
            s.fValues              = lsFun{3};
            s.fun                  = LagrangianFunction(s);
            ls1                    = DesignVariable.create(s);
            obj.designVariable{3} = ls1;
        end

        function plotMesh(obj)
            figure(1); clf;
            obj.meshSwan.plot;
            % figure(1); clf;
            % m = obj.mesh;
            % pdeplot(m.p,m.e,m.t,'xydata',m.t(4,:),'xystyle','flat','colormap','gray',...
            %               'xygrid','off','colorbar','off','title','Geometry'); 
            % axis image; axis off;
            % 
            % pdemesh(m.p,m.e,m.t); axis image; axis off;
                
        end

        function computeFreeNodes(obj)
            m         = obj.mesh;
            obj.phold = []; % freeze nodes
            
            for i = 1:size(m.ghold,1)
                obj.phold = cat(1,obj.phold,unique(m.t(1:3,m.t(4,:)==m.ghold(i))));
            end
            
            obj.pfree = setdiff(1:size(m.p,2),obj.phold); % free nodes
        end

        function solveInitialElasticProblem(obj)
            m                  = obj.mesh;
            psi_hold_all       = ones(length(m.p),obj.nMat-1);
            psi_hold_all(:,1)  = -1;
            
            obj.psi            = psi_hold_all;
            
            C = computeElasticTensor(obj);
            obj.tensor = C;

            [obj.displacements,obj.forces] = solveFEM(obj);
            obj.volume                     = computeVolume(obj);
        end

        function createOptimizationParameters(obj)
            F                       = obj.forces;
            U                       = obj.displacements;
            obj.optParams.energy0   = 0.5 * dot(F,U); % initial energy 
            max_vol                 = 2;
            obj.optParams.max_vol   = max_vol; 
            obj.optParams.voltarget = max_vol.*obj.params.volfrac; 
            obj.optParams.volstop   = obj.params.voleps.*max_vol;
            
            if obj.params.penalization == 1
                obj.optParams.volstop = max_vol;
            end
        end

        function createMassMatrix(obj)
            p = obj.mesh.p;
            t = obj.mesh.t;

            [~,obj.unitM,~] = assema(p,t,0,1,0); % mass matrix of unity density
        end

        % function createMassMatrixSwan(obj)
        %     s.test  = LagrangianFunction.create(obj.meshSwan,1,'P1');
        %     s.trial = LagrangianFunction.create(obj.meshSwan,1,'P1');
        %     s.mesh  = obj.meshSwan;
        %     s.type  = 'MassMatrix';
        %     LHS = LHSintegrator.create(s);
        %     obj.unitM = LHS.compute;     
        % end

        function updateDesignVariable(obj)
            
            levelSet = obj.psi; 
            levelSet = levelSet./normL2( obj.unitM,levelSet ); % level-set function nomalization  
            obj.psi  = levelSet;
            lsFun{1} = levelSet(:,1);
            lsFun{2} = levelSet(:,2);
            lsFun{3} = levelSet(:,3);
            
            obj.designVariable{1}.update(lsFun{1});
            obj.designVariable{2}.update(lsFun{2});
            obj.designVariable{3}.update(lsFun{3});
            
        end

        function computeThetaAndCosinus(obj)
            for i = 1 : (obj.nMat-1)
                cos(i) = obj.psi(:,i)'*obj.unitM*obj.DT(:,i);
            end
            cos            = max(min(sum(cos),1.0),-1.0);
            obj.cosAngle   = cos;

            theta          = max(real(acos(cos)),1.0e-4);
            obj.thetaAngle = theta;
        end

        function updatePlotAndDisplay(obj)
            disp(['iter    = ', num2str(obj.iter)]);
            disp(['volume  = ', num2str(obj.volume(1:end-1)),' => ', sprintf(repmat(' %3.2f%%',1,length(obj.volume)-1),obj.volume(1:end-1)*100/obj.optParams.max_vol)]);
            disp(['sf      = ', num2str(obj.shapeFunc)]);
            disp(['k       = ', num2str(obj.k)]);
            disp(['theta   = ', num2str(obj.thetaAngle*180/pi)]);
            disp(['penalty = ', num2str(obj.params.penalty)]);
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
                    
            facecolor = [186 212 244]/255;
            linecolor = [0 62 102]/255;

            p          = obj.mesh.p;
            t          = obj.mesh.t;
            fi         = obj.fiFunction;

            figure(3)       
            multimat_plot( p,t,fi );
            drawnow
        end

        function createMonitoring(obj)
            titles  = {'Vol Constraint 1'; 'Vol Constraint 2'; 'Vol Constraint 3'; 'Vol Constraint void'; 'Line Search';'Compliance';'Theta'};
            chartTypes = {'plot'; 'plot'; 'plot'; 'plot'; 'bar'; 'plot'; 'plot'};
            
            s.shallDisplay = true;
            s.maxNColumns = 7;
            s.titles = titles;
            s.chartTypes = chartTypes;

            obj.monitoring = Monitoring(s);

        end

        function obj = updateMonitoring(obj)
            for i=1:3
                vol_constraint(i) = (obj.volume(i)/obj.optParams.voltarget(i))-1;
            end

            voltarget_void = 0.8; % Beam
            %voltarget_void = 2.4; % Bridge

            volume_void = 2-(obj.volume(1)+obj.volume(2)+obj.volume(3)); % Beam
            %volume_void = 6-(obj.volume(1)+obj.volume(2)+obj.volume(3)); %Bridge
            void_constraint = (volume_void/voltarget_void)-1;

            data = [vol_constraint(1); vol_constraint(2); vol_constraint(3); void_constraint];
            data = [data;obj.k];
            data = [data;obj.Epot];
            data = [data;obj.thetaAngle];
            
            obj.monitoring.update(obj.iter,data);
        end

        function increaseLagrangianMultiplier(obj)

            difvol = obj.volume(1:end-1)-obj.optParams.voltarget;
            ic = (abs(difvol) > obj.optParams.volstop); %index control
           
            if any(ic==1) 
                if obj.params.penalization == 2 % increase the penalization parameter
                    obj.params.penalty = 2.0 * obj.params.penalty;
                    obj.k = 1;
                elseif obj.params.penalization == 3 % increase the lagrangian multiplier
                    coef = obj.volume(1:end-1)./obj.optParams.voltarget; coef = coef(ic); tau = obj.params.auglag;
                    penalty = obj.params.penalty(ic);
                    penalty = penalty + (tau(ic)./obj.params.auglag(ic)) .* (max(0,penalty - obj.params.auglag(ic).*(1.0-coef))-penalty);
                    obj.params.penalty(ic) = penalty;
                    obj.k = 1;
                end
            end

        end

        function computeComplianceAndGradient(obj)
            s.designVariable = obj.designVariable; 
            s.forces = obj.forces;
            s.displacements = obj.displacements;
            s.energy0 = obj.optParams.energy0;
            s.nMat = obj.nMat;
            s.dC = obj.dC;
            s.strain = obj.strain;
            s.tgamma = obj.tgamma;

            complianceDT = ComplianceFunctionalComputer(s);
            obj.compliance = complianceDT.J;
            obj.dJCompliance = complianceDT.dJ;
        end
        
        function TD = computeVolumeConstraintinDT(obj)
            pen = obj.params.penalty;
            maxVol = obj.optParams.max_vol;
            volt = maxVol*obj.params.volfrac;
            coef = obj.volume(1:end-1) ./ volt;
            augmentedLagr = obj.params.auglag;
            nmat = obj.nMat;
            TD = obj.dJCompliance;

            for i = 1:obj.nMat
                for j = 1:obj.nMat
                    if obj.params.penalization == 1
                        if i==j
                            TD{i,j} = 0;
                        elseif j == nmat
                            TD{i,j} = TD{i,j}- pen(i)/maxVol;
                        elseif i == nmat
                            TD{i,j} = TD{i,j} + pen(j)/maxVol;
                        else
                            TD{i,j} = TD{i,j} + (pen(j) - pen(i))/maxVol;
                        end
                    elseif obj.params.penalization == 3
                        if i==j
                            TD{i,j} = 0;
                        elseif j == nmat
                            TD{i,j} = TD{i,j} - ( (pen(i) + augmentedLagr(i)*(coef(i)-1))/volt(i) );
                        elseif i == nmat
                            TD{i,j} = TD{i,j} + ( (pen(j) + augmentedLagr(j)*(coef(j)-1))/volt(j) );
                        else
                            TD{i,j} = TD{i,j} + ( (pen(j) + augmentedLagr(j)*(coef(j)-1))/volt(j) )...
                                                      - ( (pen(i) + augmentedLagr(i)*(coef(i)-1))/volt(i) );
                        end
                    end
                end
            end
 
        end

        function dt = smoothTopologicalDerivative(obj,TD)

            s.psi = obj.psi;
            s.designVariable = obj.designVariable;
            s.m = obj.meshSwan;
            s.p = obj.mesh.p;
            s.t = obj.mesh.t;
  
            charfun = CharacteristicFunctionComputer(s); 
            [~,tfi] = charfun.computeFiandTfi();

            t = obj.mesh.t;
            p = obj.mesh.p;
            [tXi2,~] = integ_exact(t,p,obj.psi(:,2)); chi2 = (1 - tXi2); %- Mixed formulation method
            [tXi3,~] = integ_exact(t,p,obj.psi(:,3)); chi3 = (1 - tXi3); %- Mixed formulation method
            %     fi = (pdeintrp(p,t,fi)).'; % interpolation at gauss point - P1 projection method
            %     chi2 = (pdeintrp(p,t,(psi(:,2)<0))).'; %- P1 projection method
    
            dt = [];
            dt(1,:) = - tfi(1,:).*TD{1,end} - tfi(2,:).*TD{2,end} - tfi(3,:).*TD{3,end} ...
                + tfi(4,:).*( (1-chi2).*TD{end,1} + (1-chi3).*chi2.*TD{end,2} + chi2.*chi3.*TD{end,3} );
            dt(2,:) = - tfi(2,:).*TD{2,1} - tfi(3,:).*TD{3,1} + tfi(1,:).*( (1-chi3).*TD{1,2} + chi3.*TD{1,3} );
            dt(3,:) = tfi(2,:).*TD{2,3} - tfi(3,:).*TD{3,2};

            dt = pdeprtni(p,t,dt);
        end
        
        


    end

        
end