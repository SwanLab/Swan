classdef AcademicMonitoring < handle
    
    properties (Access = public)
        
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
                case 'AugmentedLagrangian'
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
            hold on
            subplot(1,3,1)
            plot(x(1),x(2),'r')
            xlabel('x_1')
            ylabel('x_2')
            title('Design variables')
            subplot(1,3,2)
            plot(iter,c,'b')
            xlabel('Iteration')
            ylabel('Objective function J(x)')
            title('Objective function evolution')
            subplot(1,3,3)
            plot(iter,cnst(1),'b',iter,cnst(2),'g')
            xlabel('Iteration')
            ylabel('Constraint violation')
            title('Constraints violation evolution')
            legend('Constraint 1, Constraint 2')
            drawnow
        end
        
        function plotNullSpace(obj,cParams)
            disp(obj.designVariable.value);
            disp(obj.cost.value);
        end
        
        function plotAugmentedLagrangian(obj,x,cParams)
            
        end
        
        function plotBisection(obj,x,cParams)
            
        end
        
        function plotIPOPT(obj,x,cParams)
            
        end
        
        function plotMMA(obj,x,cParams)
            
        end
                
    end
    
end