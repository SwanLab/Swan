classdef OrientationAsDesignVariable < handle
    
    properties (GetAccess = public, SetAccess = private)
        cost
    end
    
    properties (GetAccess = private, SetAccess = public)
        nStressFigure
    end
    
    properties (Access = private)
        nCostFigure
        orientationUpdaterType
        phyP
        hC
        designVariable
        microVariables
        costI
        cost0
        orientationUpdater
        plotting
        iter
        niter
    end
    
    methods (Access = public)
        
        function obj = OrientationAsDesignVariable(cParams)
            obj.init(cParams);
            obj.computeInitialCost();
        end
        
        function compute(obj)
            for i = 2:obj.niter
                obj.iter = i;
                obj.updateAlpha();
                obj.updateCost();
                obj.updateStress();
                obj.updateCost();
            end
        end
        
        
        function f = plotStress(obj)
            xd = obj.phyP.mesh.coord(:,1);
            yd = obj.phyP.mesh.coord(:,2);
            connec = obj.phyP.mesh.connec;
            [nelem,nnode] = size(connec);
            xn = zeros(nelem,nnode);
            yn = zeros(nelem,nnode);
            for inode = 1:nnode
                nodes = connec(:,inode);
                xn(:,inode) =  xd(nodes);
                yn(:,inode) =  yd(nodes);
            end
            xp = mean(xn,2);
            yp = mean(yn,2);
            alpha = obj.designVariable.alpha;
            f = figure(obj.nStressFigure);
            quiver(xp,yp,alpha(1,:)',alpha(2,:)') ;
            drawnow
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.plotting = cParams.plotting;
            obj.orientationUpdaterType = cParams.orientationUpdaterType;
            obj.createPhysicalProblem();
            obj.createMicroVariables();
            obj.createInitialOrientation();
            obj.createHomogenizerVarComputer();
            obj.createOrientationUpdater();
            obj.nStressFigure = 200;
            obj.nCostFigure = 100;
            obj.niter = 10;
        end
        
        function updateAlpha(obj)
            cParams.pD = obj.phyP.variables.principalDirections;
            cParams.pS = obj.phyP.variables.principalStress;
            obj.orientationUpdater.compute(cParams);
            alpha = obj.orientationUpdater.alpha;
            obj.designVariable.alpha = alpha;
            obj.updateHomogenizedMaterialProperties();
        end
        
        function updateCost(obj)
            obj.computeCost();
            obj.plotAndPrint();
        end
        
        function updateStress(obj)
            obj.updateHomogenizedMaterialProperties();
            obj.solvePDE();
        end
        
        function plotAndPrint(obj)
            if obj.plotting
                obj.plotCost();
                obj.plotStress();
                disp(obj.cost(obj.iter));
            end
        end
        
        function createHomogenizerVarComputer(obj)
            s.type = 'ByVademecum';
            s.fileName = 'SmoothRectangle';
            s.designVariable = obj.designVariable;
            s.nelem = size(obj.phyP.mesh.connec,1);
            s.designVariable = obj.designVariable;
            obj.hC = HomogenizedVarComputer.create(s);
        end
        
        function createMicroVariables(obj)
            nelem = size(obj.phyP.mesh.connec,1);
            m1 = 0.8*ones(nelem,1);
            m2 = 0.8*ones(nelem,1);
            obj.microVariables = [m1 m2];
        end
        
        function createOrientationUpdater(obj)
            cParams.type = obj.orientationUpdaterType;
            obj.orientationUpdater = OrientationUpdater.create(cParams);
        end
        
        function createInitialOrientation(obj)
            nelem = size(obj.phyP.mesh.connec,1);
            s =  SettingsDesignVariable();
            s.type = 'MicroParams';
            s.mesh = obj.phyP.mesh;
            s.scalarProductSettings.epsilon = 1e-3;
            obj.designVariable = DesignVariable.create(s);
            obj.designVariable.alpha = zeros(2,nelem);
            obj.designVariable.alpha(1,:) = 1;
            obj.designVariable.alpha(2,:) = 0;
        end
        
        function computeInitialCost(obj)
            obj.iter = 1;
            obj.updateHomogenizedMaterialProperties();
            obj.solvePDE();
            obj.cost0 = 1;
            obj.computeCost();
            obj.cost0 = obj.cost;
            obj.cost(obj.iter) = 1;
            obj.plotAndPrint();
        end
        
        function computeCost(obj)
            c = obj.phyP.computeCompliance();
            obj.cost(obj.iter) = c/obj.cost0;
        end
        
        function solvePDE(obj)
            obj.phyP.computeVariables();
        end
        
        function createPhysicalProblem(obj)
            testName = 'BridgeArchSmall';
            obj.phyP = FEM.create(testName);
            obj.phyP.preProcess();
        end
        
        function updateHomogenizedMaterialProperties(obj)
            obj.hC.computeCtensor(obj.microVariables);
            obj.phyP.computeStiffnessMatrix();
            obj.phyP.setC(obj.hC.C);
        end
        
        function plotCost(obj)
            figure(obj.nCostFigure);
            plot(obj.cost);
            drawnow
        end
        
    end
    
end