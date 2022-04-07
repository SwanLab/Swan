classdef IterativeProcessComputer < handle 

     properties (Access = public)
         designVariable
         dualVariable
         optimizerType
     end

     properties (Access = private)
         optimizer
         cost
         constraint
         stiffnessMatComputer
         bendingMatComputer
         eigenModes
     end
    
    properties (Access = private)
        nElem
        nConstraints
        length
        nValues
        youngModulus
        inertiaMoment
        minThick
        maxThick
        nIter
        maxIter
        e
        E1
        E2
        freeNodes
    end

    properties (Access = private)
        mmaVarComputer
        xMin
        xMax
        xOld1
        xOld2
        lOW
        uPP
        a0Val
        aMMA
        dVal
        cVal
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
            obj.computeBoundaryConditions();            
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
            obj.nElem          = cParams.nElem;
            obj.nConstraints   = cParams.nConstraints;
            obj.length         = cParams.length;
            obj.youngModulus   = cParams.youngModulus;
            obj.inertiaMoment  = cParams.inertiaMoment;
            obj.minThick       = cParams.minThick;
            obj.maxThick       = cParams.maxThick;
            obj.nValues        = cParams.nValues;
            obj.nIter          =  0;
            obj.mmaVarComputer = cParams.mmaVarComputer;
            obj.maxIter        = cParams.maxIter;
            obj.optimizerType  = cParams.optimizerType;
            obj.designVariable = cParams.designVariable;
         end



         function obj = computeIterativeProcess(obj)
% Refactor Constraint
% Construct Optimizer;
%              s.designVar = obj.designVariable;
%              s.type     = obj.optimizerType;
%              s.constraintCase = 'INEQUALIY';
%              s.cost = obj.cost;
%              s.constraint = obj.constraint;
%              s.dualVariable.value = zeros(1,obj.nConstraints); 
%              s.maxIter = obj.maxIter;
%              s.incrementalScheme.iStep = 1;
%              s.incrementalScheme.nSteps = 1;
%              s.targetParameters.optimality_tol = 0.0005;
%              s.historyPrinterSettings = [];
%              s.uncOptimizerSettings.ub = 10;
%              s.uncOptimizerSettings.lb = 0.25;
%              s.historyPrinterSettings.shallPrint = false;
%              s.historyPrinterSettings.fileName = 'OptimalBuckling';
%              s.optimizerNames.type = obj.optimizerType;
% 
% 
%             sm.showOptParams         = true;
%             sm.refreshInterval       = 1;
%             sm.problemID             = [];
%             sm.costFuncNames         = [];
%             sm.costWeights           = [];
%             sm.constraintFuncs       = [];
% 
%             sm.shallDisplayDesignVar = false;
%             sm.shallShowBoundaryConditions = [];
%             sm.boundaryConditions = [];
%             sm.designVariable = obj.designVariable;
%             sm.optimizerNames.type = obj.optimizerType;
%             sm.dim = [];
%             sm.scale = [];
%             sm.mesh = [];
% 
%              s.monitoringDockerSettings = sm;
% 
%              sp.shallPrint = false;
%              s.postProcessSettings = sp;
% 
%              obj.optimizer = Optimizer.create(s);    
%              obj.optimizer.solveProblem();

             obj.change = 1;
             obj.hasFinished = 0;
             while ~obj.hasFinished
                obj.increaseIter();
                obj.updateStatus();
                obj.computeNewDesign();
                obj.displayIteration()
                obj.plotFigures();
            end

         end

        function increaseIter(obj)
            obj.nIter = obj.nIter+1;
        end

        function updateStatus(obj)
            obj.hasFinished = (obj.change <= 0.0005) || (obj.nIter >= obj.maxIter);
        end

        function computeBoundaryConditions(obj)
            N = obj.nElem;
            fixnodes = union([1,2], [2*N+1,2*N+2]);
            nodes = 1:2*N+2;
            free  = setdiff(nodes,fixnodes);
            obj.freeNodes = free;
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
            sF.nElem          = obj.nElem;
            sF.type = 'firstEignValue_functional';            
            sC.weights = 1;
            sC.nShapeFuncs = 1;
            sC.designVar = obj.designVariable;         
            sC.shapeFuncSettings{1} = sF;
            obj.cost = Cost(sC);
        end

        function createEigModes(obj)
            s.freeNodes  = obj.freeNodes;
            s.nElem      = obj.nElem;
            s.length     = obj.length;
            s.stiffnessMatComputer = obj.createStiffnessMatrix();
            s.bendingMatComputer   = obj.createBendingMatrix();
            s.designVariable       = obj.designVariable;
            obj.eigenModes = EigModes(s);
        end

        function K = createStiffnessMatrix(obj)
            s.nElem          = obj.nElem;
            s.length         = obj.length;
            s.youngModulus   = obj.youngModulus;
            s.inertiaMoment  = obj.inertiaMoment;
            s.freeNodes      = obj.freeNodes;
            K = StiffnessMatrixComputer(s);
        end

        function B = createBendingMatrix(obj)
            s.nElem          = obj.nElem;
            s.length         = obj.length;
            s.youngModulus   = obj.youngModulus;
            s.inertiaMoment  = obj.inertiaMoment;
            s.freeNodes      = obj.freeNodes;
            s.designVariable = obj.designVariable;
            B = BendingMatrixComputer(s);
        end



        function createConstraint(obj)
            sF1.nElem          = obj.nElem;
            sF1.eigModes       = obj.eigenModes;
            sF1.eigNum         = 1;
            sF1.type = 'doubleEig';  

            sF2.nElem          = obj.nElem;
            sF2.eigModes       = obj.eigenModes;  
            sF2.eigNum         = 2;
            sF2.type = 'doubleEig';  

            sF3.nElem          = obj.nElem;
            sF3.type = 'volumeColumn';              


            sC.nShapeFuncs = 3;
            sC.designVar = obj.designVariable;   
            sC.dualVariable  = [];
            sC.shapeFuncSettings{1} = sF1;
            sC.shapeFuncSettings{2} = sF2;
            sC.shapeFuncSettings{3} = sF3;
            obj.constraint = Constraint(sC);
        end

        function displayIteration(obj)
            x = obj.designVariable.value;
            N = obj.nElem;
            f0val = obj.cost.value;
            iter  = obj.nIter;
            disp([' It.: ' sprintf('%4i',iter) ' Obj.: ' sprintf('%10.4f',f0val) ...
                ' Vol.: ' sprintf('%6.3f',  (1/N)*(sum(x)-x(N+1))  ) ...
                ' ch.: ' sprintf('%6.3f',''  )])
            %(obj.constraint.D(2,2)-obj.constraint.D(1,1))
        end

        function plotFigures(obj)
            N = obj.nElem;
            iter = obj.nIter;
            x = obj.designVariable.value;
            obj.costHistory(iter) = obj.cost.value;
            obj.vol(iter) = (1/N)*sum(x(1:N));
            obj.plot(iter)
            figure(3)
            plot(obj.costHistory)
            grid on
            grid minor
            xlabel('Number of Iteration','Interpreter', 'latex','fontsize',18,'fontweight','b');
            ylabel('Cost','Interpreter', 'latex','fontsize',18,'fontweight','b');
            figure(4)
            plot(obj.vol)
            grid on
            grid minor
            xlabel('Number of Iteration','Interpreter', 'latex','fontsize',18,'fontweight','b');
            ylabel('Volume','Interpreter', 'latex','fontsize',18,'fontweight','b');
        end

        function plot(obj,iter)
            x = obj.designVariable.value;
            obj.eigenModes.plot(x,iter);
        end                            
           
    end

end