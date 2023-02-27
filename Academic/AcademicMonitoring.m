classdef AcademicMonitoring < handle
    
    properties (Access = private)
        xVec1 = []
        xVec2 = []
        iterVec = []
        cVec = []
        cnstVec1 = []
        cnstVec2 = []
    end
    
    properties (Access = private)
        optimizerName
        incrementalScheme
        cost
        constraint
        designVariable
        dualVariable
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = AcademicMonitoring(cParams)
            obj.init(cParams);
        end
        
        function compute(obj,cParams)
            switch obj.optimizerName
                case 'fmincon'
                    obj.plotFmincon(cParams);
                case 'NullSpace'
                    obj.plotNullSpace(cParams);
                case 'AlternatingPrimalDual'
                    obj.plotAugmentedLagrangian(cParams);
                case 'DualNestedInPrimal'
                    obj.plotBisection(cParams);
                case 'IPOPT'
                    obj.plotIPOPT(cParams);
                case 'MMA'
                    obj.plotMMA(cParams);
                otherwise
                    error('Optimizer not implemented')
            end
        end

        function create(obj,cParams)

        end
        
    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.optimizerName     = cParams.type;
            obj.cost              = cParams.cost;
            obj.constraint        = cParams.constraint;
            obj.designVariable    = cParams.designVar;
            obj.dualVariable      = cParams.dualVariable;
            obj.incrementalScheme = cParams.incrementalScheme;
        end
        
        function plotFmincon(obj,cParams)
            iter = cParams.nIter;
            x    = obj.designVariable.value;
            c    = obj.cost.value;
            cnst = obj.constraint.value;
            obj.xVec1 = [obj.xVec1;x(1)];
            obj.xVec2 = [obj.xVec2;x(2)];
            obj.iterVec = [obj.iterVec;iter];
            obj.cVec = [obj.cVec;c];
            obj.cnstVec1 = [obj.cnstVec1;cnst(1)];
            obj.cnstVec2 = [obj.cnstVec2;cnst(2)];
            subplot(1,3,1)
            hold on
            plot(obj.xVec1,obj.xVec2,'r')
            hold off
            xlabel('x_1')
            ylabel('x_2')
            title('Design variables')
            subplot(1,3,2)
            hold on
            plot(obj.iterVec,obj.cVec,'b')
            hold off
            xlabel('Iteration')
            ylabel('Objective function J(x)')
            title('Objective function evolution')
            subplot(1,3,3)
            hold on
            plot(obj.iterVec,obj.cnstVec1,'b',obj.iterVec,obj.cnstVec2,'g')
            hold off
            xlabel('Iteration')
            ylabel('Constraint violation')
            title('Constraints violation evolution')
            legend('Constraint 1', 'Constraint 2')
            drawnow
        end
        
        function plotNullSpace(obj,cParams)
            disp(obj.designVariable.value);
            disp(obj.cost.value);
        end
        
        function plotAugmentedLagrangian(obj,x,cParams)
            disp(obj.designVariable.value);
            disp(obj.cost.value);   
        end
        
        function plotBisection(obj,x,cParams)
            
        end
        
        function plotIPOPT(obj,x,cParams)
            
        end
        
        function plotMMA(obj,x,cParams)
            
        end
                
    end
    
end