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
        end

        function compute(obj)
            obj.computeIterativeProcess();
        end
    end

    methods (Access = private)

         function init(obj,cParams)
            obj.freeNodes      = cParams.data.freeNodes;
            obj.nConstraints   = cParams.nConstraints;
            obj.mesh           = cParams.data.mesh;
            obj.dim            = cParams.data.dim;
            obj.youngModulus   = cParams.youngModulus;
            obj.inertiaMoment  = cParams.inertiaMoment;
            obj.nValues        = cParams.nValues;
            obj.nIter          =  0;
            obj.mmaVarComputer = cParams.mmaVarComputer;
            obj.maxIter        = cParams.maxIter;
            obj.optimizerType  = cParams.optimizerType;
            obj.designVariable = cParams.bProblem.designVariable;
            obj.cost           = cParams.bProblem.cost;
            obj.constraint     = cParams.bProblem.constraint;
         end

         function obj = computeIterativeProcess(obj)
% -------
            s = SettingsOptimizer();
            s.optimizerNames.type = obj.optimizerType; %'MMA';% 'AlternatingPrimalDual';'NullSpace';'MMA';'AlternatingPrimalDual';'NullSpace';'AlternatingPrimalDual';'MMA';'AlternatingPrimalDual';%MMA';%'IPOPT';%fmincon';%'MMA';%'fmincon';'MMA';

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
            sD.nConstraints = 3;
            s.dualVariable     = DualVariable(sD);              
            s.uncOptimizerSettings.ub = obj.mmaVarComputer.maxThick;
            s.uncOptimizerSettings.lb = obj.mmaVarComputer.minThick;        
            s.outputFunction.type        = 'Topology';
            s.outputFunction.iterDisplay = 'none';
            s.type = obj.optimizerType;%'MMA';% 'AlternatingPrimalDual';%'fmincon';'AlternatingPrimalDual';'NullSpace';'AlternatingPrimalDual';'MMA';'AlternatingPrimalDual';'MMA'; % IPOPT';%'fmincon';'MMA';%'fmincon';%'MMA';
            s.outputFunction.monitoring  = MonitoringManager(s);                  
            s.maxIter           = obj.maxIter;
            s.constraintCase = {'INEQUALITY','INEQUALITY','INEQUALITY'};
            %s.primalUpdater = 'PROJECTED GRADIENT';

            obj.optimizer = Optimizer.create(s);    
            obj.optimizer.solveProblem();
         end

        function plot(obj,iter)
            A = obj.designVariable.getColumnArea;
            obj.eigenModes.plot(A,iter);
        end 


                                   
           
    end

end
