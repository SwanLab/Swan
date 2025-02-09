classdef IterativeProcessComputer < handle 

     properties (Access = public)
         designVariable
         dualVariable
         optimizerType
     end

     properties (Access = private) % Classes
         mesh
         dim
         optimizer
         cost
         constraint
         eigenModes
     end
    
    properties (Access = private)
        nConstraints
        nValues
        youngModulus
        inertiaMoment
        nIter
        maxIter
        freeNodes
    end

    properties (Access = private)
        mmaVarComputer
        hasFinished
        change
        costHistory
        xMMA
        xVal
        vol
    end

    methods (Access = public)

        function obj = IterativeProcessComputer(cParams)
            obj.init(cParams);
            obj.createCost();
            obj.createEigModes();
            obj.createConstraint();
        end

        function compute(obj)
            obj.computeIterativeProcess();
        end
    end

    methods (Access = private)

         function init(obj,cParams)
            obj.freeNodes      = cParams.freeNodes;
            obj.nConstraints   = cParams.nConstraints;
            obj.mesh           = cParams.mesh;
            obj.dim            = cParams.dim;
            obj.youngModulus   = cParams.youngModulus;
            obj.inertiaMoment  = cParams.inertiaMoment;
            obj.nValues        = cParams.nValues;
            obj.nIter          =  0;
            obj.mmaVarComputer = cParams.mmaVarComputer;
            obj.maxIter        = cParams.maxIter;
            obj.optimizerType  = cParams.optimizerType;
            obj.designVariable = cParams.designVariable;
         end

         function obj = computeIterativeProcess(obj)
% -------
            s = SettingsOptimizer();
            s.optimizerNames.type = 'NullSpace';'AlternatingPrimalDual';'NullSpace';'AlternatingPrimalDual';'MMA';'AlternatingPrimalDual';%MMA';%'IPOPT';%fmincon';%'MMA';%'fmincon';'MMA';
            s.optimizerNames.primal = 'PROJECTED GRADIENT';
            s.uncOptimizerSettings.scalarProductSettings = obj.designVariable.scalarProduct;
            s.uncOptimizerSettings.designVariable   = obj.designVariable;
            s.monitoringDockerSettings.mesh = obj.mesh;
            s.monitoringDockerSettings.optimizerNames = s.optimizerNames;
            s.monitoringDockerSettings.refreshInterval = 1;
            s.designVar         = obj.designVariable;
            s.targetParameters.optimality_tol  = 0.0005; %obj.incrementalScheme.targetParams;
            s.targetParameters.constr_tol = 0.0005;
            s.cost              = obj.cost;
            s.constraint        = obj.constraint;
            s.incrementalScheme.iStep  = 1;%obj.incrementalScheme;
            s.incrementalScheme.nSteps = 1;
            s.dualVariable      = obj.dualVariable;              
            s.uncOptimizerSettings.ub = 10;
            s.uncOptimizerSettings.lb = 0.25;        
            s.outputFunction.type        = 'Topology';
            s.outputFunction.iterDisplay = 'none';
            s.type = 'NullSpace';'AlternatingPrimalDual';'NullSpace';'AlternatingPrimalDual';'MMA';'AlternatingPrimalDual';'MMA';%IPOPT';%'fmincon';'MMA';%'fmincon';%'MMA';
            s.outputFunction.monitoring  = MonitoringManager(s);                  
            s.maxIter           = 1000;
            s.constraintCase = 'INEQUALITY';
            %s.primalUpdater = 'PROJECTED GRADIENT';

            obj.optimizer = Optimizer.create(s);    
            obj.optimizer.solveProblem();
% -------
%              obj.change = 1;
%              obj.hasFinished = 0;
%              while ~obj.hasFinished
%                 obj.increaseIter();
%                 obj.updateStatus();
%                 obj.computeNewDesign();
%                 obj.displayIteration()
%                 obj.plotFigures();
%             end
%
         end

         function increaseIter(obj)
             obj.nIter = obj.nIter+1;
         end

         function updateStatus(obj)
             obj.hasFinished = (obj.change <= 0.0005) || (obj.nIter >= obj.maxIter);
        end

        function computeNewDesign(obj)
            iter = obj.nIter;
            x = obj.designVariable.value;
            xval = x;

            obj.constraint.computeFunctionAndGradient(); 
            fval = obj.constraint.value;
            dfdx = obj.constraint.gradient';
            dfdx2 = 0;            
            
            obj.cost.computeFunctionAndGradient(); 
            f0val = obj.cost.value;
            df0dx = obj.cost.gradient;
            df0dx2 = 0;
            
            nV = obj.nValues;
            m = obj.nConstraints;
            xmma = obj.mmaVarComputer.compute(nV,m,iter,xval,f0val,df0dx,df0dx2,fval,dfdx,dfdx2);

            obj.xVal = xval;
            obj.xMMA = xmma;
  
            obj.designVariable.update(obj.xMMA);
            obj.change = max(abs(obj.designVariable.value-xval));            
        end

        function createCost(obj)
            sF.type = 'firstEignValue_functional';            
            sC.weights = 1;
            sC.nShapeFuncs = 1;
            sC.designVar = obj.designVariable;         
            sC.shapeFuncSettings{1} = sF;
            obj.cost = Cost(sC);
        end

        function createEigModes(obj)
            s.freeNodes  = obj.freeNodes;
            s.mesh       = obj.mesh;
            s.dim        = obj.dim;
            s.stiffnessMatComputer = obj.createStiffnessMatrix();
            s.bendingMatComputer   = obj.createBendingMatrix();
            s.designVariable       = obj.designVariable;
            obj.eigenModes = EigModes(s);
        end

        function K = createStiffnessMatrix(obj)
            s.type = 'StiffnessMatrixColumn';
            s.dim = obj.dim;
            s.mesh = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.freeNodes      = obj.freeNodes;
            K = LHSintegrator.create(s);
        end

        function B = createBendingMatrix(obj)
            s.type         = 'BendingMatrix';
            s.dim          = obj.dim;
            s.mesh         = obj.mesh;
            s.globalConnec = obj.mesh.connec;
            s.inertiaMoment  = obj.inertiaMoment;
            s.youngModulus   = obj.youngModulus;
            s.designVariable = obj.designVariable;
            s.freeNodes      = obj.freeNodes;
            B = LHSintegrator.create(s);
        end

        function createConstraint(obj)
            sF1.eigModes       = obj.eigenModes;
            sF1.eigNum         = 1;
            sF1.type = 'doubleEig';  

            sF2.eigModes       = obj.eigenModes;  
            sF2.eigNum         = 2;
            sF2.type = 'doubleEig';  

            sF3.type = 'volumeColumn';    
            sF3.mesh = obj.mesh;

            sC.nShapeFuncs = 3;
            sC.designVar = obj.designVariable;   
            sC.dualVariable  = [];
            sC.shapeFuncSettings{1} = sF1;
            sC.shapeFuncSettings{2} = sF2;
            sC.shapeFuncSettings{3} = sF3;
            obj.constraint = Constraint(sC);
        end

        function displayIteration(obj)
            V = obj.designVariable.computeVolum();
            f0val = obj.cost.value;
            iter  = obj.nIter;
            disp([' It.: ' sprintf('%4i',iter) ' Obj.: ' sprintf('%10.4f',f0val) ...
                ' Vol.: ' sprintf('%6.3f',V) ...
                ' ch.: ' sprintf('%6.3f',''  )])
        end

        function plotFigures(obj)
            iter = obj.nIter;
            obj.costHistory(iter) = obj.cost.value;
            obj.vol(iter) = obj.designVariable.computeVolum();
            obj.plot(iter)
            figure(4)
            plot(obj.costHistory)
            grid on
            grid minor
            xlabel('Number of Iteration','Interpreter', 'latex','fontsize',18,'fontweight','b');
            ylabel('Cost','Interpreter', 'latex','fontsize',18,'fontweight','b');
            figure(5)
            plot(obj.vol)
            grid on
            grid minor
            xlabel('Number of Iteration','Interpreter', 'latex','fontsize',18,'fontweight','b');
            ylabel('Volume','Interpreter', 'latex','fontsize',18,'fontweight','b');
        end

        function plot(obj,iter)
            A = obj.designVariable.getColumnArea;
            obj.eigenModes.plot(A,iter);
        end                            
           
    end

end