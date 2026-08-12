classdef TutorialHomogenization < handle

    properties (Access = public)
        paramHole
        Chomog
        volFrac
        f,df,ddf
    end

    properties (Access = private)
        E
        nu
        meshType
        meshN
        holeType
        nSteps
        % damageType
        pnorm
        monitoring
        Mmass
        % lattice
        theta
        a1
        a2
        c
    end

    properties (Access = private)
        baseMesh
        masterSlave
        test
        maxParam
        
    end

    methods (Access = public)
        
        function obj = TutorialHomogenization(c1,c2,alpha,meshN)

            if nargin == 0
                % modo padrão
                obj.init();
            else
                obj.init();
                obj.meshType = 'Square';
                obj.c = [c1,c2];
                obj.theta = [0,alpha];
                obj.meshN = meshN;
            end
        
            obj.defineMesh();
        end
        function [C,vf] = runAtVolume(obj,targetVF)

            obj.computeHoleParams();
            obj.compute();
        
            [~, idx] = min(abs(obj.volFrac - targetVF));
            C = obj.Chomog(:,:,:,:,idx);
            vf = obj.volFrac(idx);
        
        end
      
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.E          = 1;
            obj.nu         = 0.3;
            obj.meshType   = 'Square';
            obj.meshN      = 40;

            obj.holeType   = 'Square';
            obj.pnorm      = 'Inf';
            % obj.damageType = 'Area';
            obj.nSteps     = 20;

            obj.monitoring = true;
        end

        function defineMesh(obj)
            switch obj.meshType
                case 'Square'
                    % s.c = [1,1];
                    % obj.theta = [0,70];
                    s.theta = obj.theta;
                    s.c = obj.c;
                    s.theta = obj.theta;
                    obj.c     = s.c;
                    obj.theta = s.theta;
                    
                    obj.a1 = obj.c(1) * [cosd(obj.theta(1)), sind(obj.theta(1))];
                    obj.a2 = obj.c(2) * [cosd(obj.theta(2)), sind(obj.theta(2))];
                    s.divUnit = obj.meshN;
                    s.filename = '';
                    MC = MeshCreator(s);
                    MC.computeMeshNodes();
                    
                    % obj.lattice = MC.lattice;
                case 'Hexagon'
                    s.c = [1.3,0.6,0.4];
                    s.theta = [0,60,120];
                    s.divUnit = obj.meshN;
                    s.filename = '';
                    MC = MeshCreator(s);
                    MC.computeMeshNodes();
                    % obj.lattice = MC.lattice;
            end
            s.coord  = MC.coord;
            s.connec = MC.connec;
            obj.baseMesh = Mesh.create(s);
            obj.masterSlave = MC.masterSlaveIndex;
            obj.test = LagrangianFunction.create(obj.baseMesh,1,'P1');
            obj.Mmass = IntegrateLHS(@(u,v) DP(v,u), obj.test, obj.test, obj.baseMesh, 'Domain');
        end

        function computeHoleParams(obj)
            obj.maxParam = 0.979*ones(size(obj.nSteps));
            nParam = length(obj.maxParam);
            obj.paramHole = cell(1,nParam);
            for i=1:nParam
                obj.paramHole{i} = linspace(1e-9,obj.maxParam(i),obj.nSteps(i));
            end
        end

        function compute(obj)
            comb = table2array(combinations(obj.paramHole{:}));
            nComb = size(comb,1);
            mat = zeros(2,2,2,2,nComb);
            volF = zeros(1,nComb);
            for i=1:nComb
                hole = comb(i,:);
                if i==1
                    hole = 1e-10*ones(size(hole));
                end
                mat(:,:,:,:,i) = obj.computeHomogenization(hole);
                volF(i)    = obj.computeVolumeFraction(hole);
            end
            obj.Chomog = obj.assembleResults(mat);
            obj.volFrac = obj.assembleResults(volF);
        end
        
        function matHomog = computeHomogenization(obj,l)
            dens = obj.createDensityLevelSet(l);
            mat  = obj.createDensityMaterial(dens);
            matHomog = obj.solveElasticMicroProblem(mat,dens);
        end

        function lsf = createDensityLevelSet(obj,l)
            ls = obj.computeLevelSet(obj.baseMesh,l);
            sUm.backgroundMesh = obj.baseMesh;
            sUm.boundaryMesh   = obj.baseMesh.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(ls);

            ls = CharacteristicFunction.create(uMesh);
            s.trial = obj.test;
            s.mesh = obj.baseMesh;
            f = FilterLump(s); 
            lsf = f.compute(ls,2);
        end

        function ls = computeLevelSet(obj,mesh,l)
           gPar.type = obj.holeType;
           gPar.pnorm = obj.pnorm;   
           % a1 = obj.lattice.a1;
           % phi = atan2(a1(2), a1(1));
           % gPar.rotation = phi;
            coord = mesh.coord;
            xmin = min(coord(:,1));
            xmax = max(coord(:,1));
            ymin = min(coord(:,2));
            ymax = max(coord(:,2));          
            center_x = (xmin + xmax)/2;
            center_y = (ymin + ymax)/2;
            gPar.xCoorCenter = center_x;
            gPar.yCoorCenter = center_y;
          
            switch obj.meshType
                case 'Square'
                    % % gPar.xCoorCenter = 0.5;
                    % % gPar.yCoorCenter = 0.5;
                    % % a1 = obj.lattice.a1;
                    % % a2 = obj.lattice.a2;
                    % 
                    % % alpha = 0.5;
                    % % beta  = 0.5;
                    % 
                    % % center = alpha*a1 + beta*a2;
                    % 
                    % gPar.xCoorCenter = 0.5;
                    % gPar.yCoorCenter = 0.5;
                    % if size(coord,1) >= 4
                    %     v1 = coord(1,:);
                    %     v2 = coord(2,:);
                    %     v3 = coord(3,:);          
                    %     side1 = v2 - v1;
                    %     phi = atan2(side1(2), side1(1));
                    % else
                    %     phi = 0;
                    % end
                     phi = atan2(obj.a1(2), obj.a1(1));
                     gPar.a1 = obj.a1;
                     gPar.a2 = obj.a2;
                case 'Hexagon'
                    
                    if size(coord,1) >= 6
                        v1 = coord(1,:);
                        v2 = coord(2,:);
                        side1 = v2 - v1;
                        phi = atan2(side1(2), side1(1));
                    else
                        phi = 0;
                    end           
            end
            gPar.rotation = phi;
            
            switch obj.holeType
                case 'Circle'
                    gPar.radius = l/2;
                case 'Square'
                    gPar.length = l;
                    
                case 'SmoothRectangle'
                   % a1 = obj.lattice.a1(:);
                   % 
                   % 
                   %  Lc = norm(a1);  
                   % 
                   % 
                   %  ratio = 1.2 / 0.6;
                   % 
                   % 
                   %  sx = l(1) * Lc;
                   %  sy = sx / ratio;
                   % 
                   %  gPar.xSide = sx;
                   %  gPar.ySide = sy;
                    % gPar.pnorm = 16;

                    sx = l;
                    sy = l/2;            
                   
                    
                    gPar.xSide = sx;
                    gPar.ySide = sy;
                    gPar.pnorm = 16;
                case 'Ellipse'
                    gPar.type = "SmoothRectangle";
                    gPar.xSide  = l(1);
                    gPar.ySide  = l(2);
                    gPar.pnorm  = 2;  
                case 'SmoothHexagon'
                    gPar.radius = l;
                    gPar.normal = [0 1; sqrt(3)/2 1/2; sqrt(3)/2 -1/2];
                case 'ReinforcedHoneycomb'
                    gPar.theta  = 1-l;                          
                    gPar.eps    = 1;                        
                    gPar.normal = [0 1; sqrt(3)/2 1/2; sqrt(3)/2 -1/2];            
                    gPar.radius = l;    
                    gPar.rotation = phi;
            end  
            g                  = GeometricalFunction(gPar);
            phiFun             = g.computeLevelSetFunction(mesh);
            lsCircle           = phiFun.fValues;
            % if l(1) <= 1e-9 && gPar.theta == 1
            %     ls = ones(size(lsCircle));
            % else
            ls = -lsCircle; 
            % end            
        end

        function mat = createDensityMaterial(obj,lsf)
            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
            s.matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(1e-6*obj.E,obj.nu,obj.baseMesh.ndim);
            s.matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(1e-6*obj.E,obj.nu);
            s.matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(obj.E,obj.nu,obj.baseMesh.ndim);
            s.matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(obj.E,obj.nu);
            mI = MaterialInterpolator.create(s);

            x{1} = lsf;
            s.mesh                 = obj.baseMesh;
            s.type                 = 'DensityBased';
            s.density              = x;
            s.materialInterpolator = mI;
            s.dim                  = '2D';
            mat = Material.create(s);
        end

        function matHomog = solveElasticMicroProblem(obj,material,dens)           
            if obj.monitoring == true
                close all
                dens.plot
                shading interp
                colormap (flipud(pink))
                drawnow
            end

            s.mesh = obj.baseMesh;
            s.material = material;
            s.scale = 'MICRO';
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions(obj.baseMesh);
            s.solverCase = DirectSolver();
            s.solverType = 'REDUCED';
            s.solverMode = 'FLUC';
            fem = ElasticProblemMicro(s);
            material.setDesignVariable({dens})
            fem.updateMaterial(material.obtainTensor())
            fem.solve();
            fem.plotMicroFields(2,50);

            totVol = obj.baseMesh.computeVolume();
            matHomog = fem.Chomog/totVol;
        end

        function bc = createBoundaryConditions(obj,mesh)
                    switch obj.meshType
                        case 'Square'
                            isBottom = @(coor) (abs(coor(:,2) - min(coor(:,2))) < 1e-12);
                            % isTop    = @(coor) (abs(coor(:,2) - max(coor(:,2))) < 1e-12);
                            % isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1))) < 1e-12);
                            % isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1))) < 1e-12);
                            % isVertex = @(coor) (isTop(coor) & isLeft(coor))    |...
                                               % (isTop(coor) & isRight(coor))   |...
                                               % (isBottom(coor) & isLeft(coor)) |...
                                               % (isBottom(coor) & isRight(coor));
                            sDir{1}.domain    = @(coor) isBottom(coor);
                            sDir{1}.direction = [1,2];
                            sDir{1}.value     = 0;
                        case 'Hexagon'
                            isBottom      = @(coor) (abs(coor(:,2) - min(coor(:,2))) < 1e-12);
                            % isTop         = @(coor) (abs(coor(:,2) - max(coor(:,2))) < 1e-12);
                            
                            coorRotY = obj.defineRotatedCoordinates(pi/3);
                            isRightBottom = @(coor) (abs(coorRotY(coor) - min(coorRotY(coor))) < 1e-12);
                            % isLeftTop     = @(coor) (abs(coorRotY(coor) - max(coorRotY(coor))) < 1e-12);
                            % coorRotY = obj.defineRotatedCoordinates(-pi/3);
                            % isLeftBottom  = @(coor) (abs(coorRotY(coor) - min(coorRotY(coor))) < 1e-12);
                            % isRightTop    = @(coor) (abs(coorRotY(coor) - max(coorRotY(coor))) < 1e-12);
                            isVertex = @(coor) (isBottom(coor) & isRightBottom(coor))  ; 
                                               % (isRightBottom(coor) & isRightTop(coor))|...
                                               % (isRightTop(coor) & isTop(coor))        |...
                                               % (isTop(coor) & isLeftTop(coor))         |...
                                               % (isLeftTop(coor) & isLeftBottom(coor))  |...
                                               % (isLeftBottom(coor) & isBottom(coor))   ;
                            sDir{1}.domain    = @(coor) isBottom(coor);
                            sDir{1}.direction = [1,2];
                            sDir{1}.value     = 0;
                    end
        
                    dirichletFun = [];
                    for i = 1:numel(sDir)
                        dir = DirichletCondition(mesh, sDir{i});
                        dirichletFun = [dirichletFun, dir];
                    end
                    s.dirichletFun = dirichletFun;
                    s.pointloadFun = [];
                    s.periodicFun  = 1; %Set to not be empty
                    s.mesh = mesh;
                    bc = BoundaryConditions(s);
                    bc.updatePeriodicConditions(obj.masterSlave);
                    
        end
        

        function coorRot = defineRotatedCoordinates(~,theta)
            x0 = 0.5; y0 = sqrt(1-0.5^2);
            coorRot = @(coor) feval(@(fun) fun(:,2),([cos(theta) sin(theta); sin(theta) cos(theta)]*(coor-[x0,y0])')');
        end


        function fracVol = computeVolumeFraction(obj,l)
            % switch obj.damageType
            %     case 'Area'
            %         switch obj.holeType
            %             case {'Circle','Square'}
            %                 phi = l^2;
            %             case {'Ellipse','Rectangle'}
            %                 phi = l(1)*l(2);
            %             case {'SmoothHexagon','Hexagon'}
            %                 perimeter = 6*l;
            %                 apothem   = sqrt(l^2 - (l/2)^2);
            %                 phi = (perimeter*apothem)/(6*sqrt(3)/2);
            %         end
            %     case 'Perimeter'
            %         switch obj.holeType
            %             case {'Circle','Square','SmoothHexagon','Hexagon'}
            %                 phi = l;
            %             case 'Ellipse'
            %                 phi = pi*(3*(l(1)+l(2))-sqrt((3*l(1)+l(2))*(l(1)+3*l(2))))/...
            %                       pi*(3*(2)-sqrt((3+1)*(1+3)));
            %             case 'Rectangle'
            %                 phi = (l(1)+l(2))/2;
            %         end
            % end

         rho = obj.createDensityLevelSet(l);    
         volDom = Integrator.compute(ConstantFunction.create(1,obj.baseMesh),obj.baseMesh,2);
         fracVol = Integrator.compute(rho,rho.mesh,2)/volDom;
         
         
       
         % rho  = lsf.fValues(:);                 
         % one  = ones(size(rho));
         % volSolid = one' * (obj.Mmass * rho);  
         % volTotal = one' * (obj.Mmass * one);   
         % rho = volSolid / volTotal;             
   
        end
        
        function [mat] = assembleResults(obj,vec)
            sizeRes = size(vec);
            mat = zeros([sizeRes(1:end-1),obj.nSteps]);
            nStepsLastParam = obj.nSteps(end);
            nCombs = sizeRes(end);
            idxVec = repmat({':'}, 1, ndims(vec));
            idxMat = repmat({':'}, 1, ndims(mat));
            for i=1:nStepsLastParam
                idxVec{end} = i:nStepsLastParam:nCombs;
                idxMat{end} = i;
                mat(idxMat{:}) = vec(idxVec{:});
            end
        end
        
        function printTensorAtVolume(obj, targetVol)
        
            % encontrar índice mais próximo
            [~, idx] = min(abs(obj.volFrac - targetVol));
        
            fprintf('\nVolume solicitado: %.4f\n', targetVol);
            fprintf('Volume encontrado:  %.4f\n\n', obj.volFrac(idx));
        
            C = obj.Chomog(:,:,:,:,idx);
        
            % ---- converter para Voigt 2D ----
            Cvoigt = zeros(3,3);
        
            Cvoigt(1,1) = C(1,1,1,1);
            Cvoigt(1,2) = C(1,1,2,2);
            Cvoigt(1,3) = C(1,1,1,2);
        
            Cvoigt(2,1) = C(2,2,1,1);
            Cvoigt(2,2) = C(2,2,2,2);
            Cvoigt(2,3) = C(2,2,1,2);
        
            Cvoigt(3,1) = C(1,2,1,1);
            Cvoigt(3,2) = C(1,2,2,2);
            Cvoigt(3,3) = C(1,2,1,2);
        
            disp('Tensor homogenizado (Voigt):')
            disp(Cvoigt)
        
            % ---- módulos equivalentes ----
            K = (Cvoigt(1,1) + Cvoigt(2,2) + 2*Cvoigt(1,2))/4;
            mu = Cvoigt(3,3);
        
            fprintf('\nBulk modulus efetivo: %.6f\n', K);
            fprintf('Shear modulus efetivo: %.6f\n', mu);
        
            % ---- Young e Poisson equivalente (2D) ----
            Eeff = mu*(3*K + mu)/(K + mu);
            nueff = (K - mu)/(K + mu);
        
            fprintf('Young efetivo: %.6f\n', Eeff);
            fprintf('Poisson efetivo: %.6f\n', nueff);
        
            % ---- heatmap ----
            figure;
            imagesc(Cvoigt)
            colorbar
            axis equal tight
            title(['Tensor homogenizado - vol = ', num2str(obj.volFrac(idx))])
        end



        function plot(obj)
            tiledlayout(1,3)
            nexttile
            hold on
            plot(obj.volFrac,squeeze(obj.Chomog(1,1,1,1,:)),'LineStyle','none','Marker','o')
            % fplot(obj.f(1,1,1,1),[0 1])
            nexttile
            hold on
            plot(obj.volFrac,squeeze(obj.Chomog(2,2,2,2,:)),'LineStyle','none','Marker','o')
            % fplot(obj.f(1,1,2,2),[0 1])
            nexttile
            hold on
            plot(obj.volFrac,squeeze(obj.Chomog(1,2,1,2,:)),'LineStyle','none','Marker','o')
            % fplot(obj.f(1,2,1,2),[0 1]) 
           
        end

        % function fitting(obj)
        %     [obj.f,obj.df,obj.ddf] = DamageHomogenizationFitter.computePolynomial(8,obj.volFrac,obj.Chomog);
        % end
        
    end

   
    
end

