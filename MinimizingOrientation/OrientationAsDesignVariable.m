classdef OrientationAsDesignVariable < handle
    
    properties (Access = private)
        phyP
        hC
        designVariable
        regDesignVariable
        filter
        cost
        costI
        cost0
    end
    
    methods (Access = public)
        
        function obj = OrientationAsDesignVariable()
            obj.createPhysicalProblem();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createHomogenizerVarComputer();
            obj.computeInitialCost();
            obj.compute();
        end
        
    end
    
    methods (Access = private)
        
        function compute(obj)
            niter = 100;
            for i = 1:niter
                obj.updateAlpha();
                %obj.updateAlpha2();
                obj.updateHomogenizedMaterialProperties();                
                
                obj.computeCost();
                obj.plotCost();
                obj.plotStress();                
                disp(obj.cost);
                obj.updateHomogenizedMaterialProperties();
                obj.solvePDE();
                obj.computeCost();
                obj.plotCost();
                obj.plotStress();                                
                disp(obj.cost);
            end
        end
        
        function plotCost(obj)
            figure(1);
            obj.costI(end+1) = obj.cost;
            plot(obj.costI);
            drawnow
        end
        
        function plotStress(obj)
            U = obj.filter.getP1fromP0(squeeze(obj.designVariable.alpha(1,:)));
            V = obj.filter.getP1fromP0(squeeze(obj.designVariable.alpha(2,:)));
            x = obj.phyP.mesh.coord(:,1);
            y = obj.phyP.mesh.coord(:,2);
            conn = obj.phyP.mesh.connec;
            for inode = 1:3
                nodes = conn(:,inode);
                xp(:,inode) =  x(nodes);
                yp(:,inode) =  y(nodes);
            end
            xp = mean(xp');
            yp = mean(yp');
            a1 = obj.designVariable.alpha(1,:);
            a2 = obj.designVariable.alpha(2,:);
            figure(200);
            quiver(xp,yp,a1,a2) ;
            drawnow
        end
        
        function createFilter(obj)
            s = SettingsFilter();
            s.quadratureOrder = 'LINEAR';
            s.designVar = obj.designVariable;
            obj.filter = FilterFactory().create(s);
            obj.filter.preProcess();
        end
        
        
        function createHomogenizerVarComputer(obj)
            s.type = 'ByVademecum';
            s.fileName = 'SmoothRectangle';
            s.designVariable = obj.designVariable;
            s.nelem = size(obj.phyP.mesh.connec,1);
            s.designVariable = obj.designVariable;
            obj.hC = HomogenizedVarComputer.create(s);
        end
        
        function createDesignVariable(obj)
            s =  SettingsDesignVariable();
            s.type = 'MicroParams';
            s.mesh = obj.phyP.mesh;
            s.scalarProductSettings.epsilon = 1e-3;
            obj.designVariable = DesignVariable.create(s);
            obj.designVariable.alpha = zeros(2,size(obj.phyP.mesh.connec,1));
            obj.designVariable.alpha(1,:) = 1;
            obj.designVariable.alpha(2,:) = 0;
        end
        
        function computeInitialCost(obj)
            obj.filterDesignVariable();
            obj.updateHomogenizedMaterialProperties();
            obj.solvePDE();                        
            obj.cost0 = 1;
            obj.computeCost();
            obj.cost0 = obj.cost;
            obj.cost = 1;
            obj.plotCost()
            obj.plotStress();
        end
        
        
        function computeCost(obj)
            c = obj.phyP.computeCompliance();
            %u = obj.phyP.variables.d_u;
            %f = obj.phyP.variables.fext;
            %c = f'*u;
            obj.cost = c/obj.cost0;
        end
        
        function solvePDE(obj)
            obj.phyP.computeVariables();
        end
        
        function createPhysicalProblem(obj)
            %testName = 'Cantilever_triangle_fine';
            %testName = 'CantileverSquare';
            %testName = 'CantileverLong';%'CantileverSquare';'Cantilever_triangle_fine';
            testName = 'BridgeArchSmall';
            obj.phyP = FEM.create(testName);
            obj.phyP.preProcess();
        end
        
        function updateAlpha(obj)
            pD = obj.phyP.variables.principalDirections;
            pS = obj.phyP.variables.principalStress;
            [~,indM] = max(abs(pS));
            for i = 1:2
                dirD = squeeze(pD(i,:,:));
                for j = 1:2
                    ind = indM == j;
                    dir(i,ind) = dirD(j,ind);
                end
            end
            obj.designVariable.alpha = squeeze(pD(:,1,:));%dir
            %obj.designVariable.alpha = dir;
        end
        
        function updateAlpha2(obj)
            p = Pedersen();
            Ctensor = obj.hC.C;
            strain = obj.phyP.variables.strain;
            nelem = size(strain,3);
            dir = zeros(2,nelem);
            for ielem = 1:nelem
                s(:,1) = squeeze(strain(1,:,ielem));
                C = Ctensor(:,:,ielem);
                p.compute(s,C);
                a = p.optimalAngle();
                dir(:,ielem) = [cos(a);sin(a)];
                disp(ielem/nelem*100)
            end
            obj.designVariable.alpha = dir;
        end
        
        
        function filterDesignVariable(obj)
            nx = length(obj.designVariable.value)/obj.designVariable.nVariables;
            x  = obj.designVariable.value;
            for ivar = 1:obj.designVariable.nVariables
                i0 = nx*(ivar-1) + 1;
                iF = nx*ivar;
                xs = x(i0:iF);
                xf(:,ivar) = obj.filter.getP0fromP1(xs);
            end
            obj.regDesignVariable = xf;
        end
        
        function updateHomogenizedMaterialProperties(obj)
            obj.hC.computeCtensor(obj.regDesignVariable);
            obj.phyP.setC(obj.hC.C);            
        end
        
    end
    
end