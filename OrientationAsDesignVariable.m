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
                obj.plotCost(i)
                obj.plotStress();
                obj.updateAlpha();
                obj.updateHomogenizedMaterialProperties();
                obj.computeCost();
                disp(obj.cost)
            end
        end
        
        function plotCost(obj,i)
            figure(1);
            obj.costI(i) = obj.cost;                
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
            figure(2);
            quiver(xp,yp,a1,a2)            
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
           s.vademecumFileName = 'SmoothRectangle';
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
            obj.cost0 = 1;
            obj.computeCost();
            obj.cost0 = obj.cost;
            obj.cost = 1;
        end
            
            
        function computeCost(obj)
            obj.solvePDE();
            u = obj.phyP.variables.d_u;
            f = obj.phyP.variables.fext;
            c = f'*u;
            obj.cost = c/obj.cost0;            
        end
        
        function solvePDE(obj)
            obj.phyP.setC(obj.hC.C);
            obj.phyP.computeVariables();
        end
        
        function createPhysicalProblem(obj)
            %testName = 'Cantilever_triangle_fine';
            %testName = 'CantileverSquare';
            %testName = 'CantileverLong';%'CantileverSquare';'Cantilever_triangle_fine';
            testName = 'BridgeArch';
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
            obj.designVariable.alpha = dir;%squeeze(pD(:,1,:));%dir
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
        end
         
    end        
    
end