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
        shallPrint
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = AcademicMonitoring(cParams)
            obj.init(cParams);
        end

        function compute(obj,cParams)
            if (obj.shallPrint)
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
                    case 'IPM'
                        obj.plotIPM(cParams);
                    otherwise
                        error('Optimizer not implemented')
                end
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
            obj.shallPrint        = cParams.shallPrint;
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
            vector1 = obj.xVec1;
            vector2 = obj.xVec2;
            iterations = obj.iterVec;
            costVec = obj.cVec;
            constrVec1 = obj.cnstVec1;
            constrVec2 = obj.cnstVec2;
%             save('DesignVariableEvolution_fmincon_IP_Test3.mat','vector1','vector2','-mat')
%              save('Cost_Const_Test3_fmincon_IP.mat','iterations',"costVec","constrVec1","constrVec2")
        end
        
        function plotNullSpace(obj,cParams)
            disp(obj.designVariable.value);
            disp(obj.cost.value);
            x    = obj.designVariable.value;
            obj.xVec1 = [obj.xVec1;x(1)];
            obj.xVec2 = [obj.xVec2;x(2)];
            vector1 = obj.xVec1;
            vector2 = obj.xVec2;
            %save('DesignVariableEvolution_NullSpace_Test2.mat','vector1','vector2','-mat')
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

        function plotIPM(obj,cParams)
            iter = cParams.nIter;
            x    = obj.designVariable.value;
            c    = obj.cost.value;
            cnst = obj.constraint.value;
            if length(x) == 2
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
                xlabel('$x_1$','Interpreter','latex')
                ylabel('$x_2$','Interpreter','latex')
                title('Design variables','Interpreter','latex')
                set(gca,'TickLabelInterpreter','latex')
                box on
                subplot(1,3,2)
                hold on
                plot(obj.iterVec,obj.cVec,'b')
                hold off
                xlabel('Iteration','Interpreter','latex')
                ylabel('Objective function $J(x)$','Interpreter','latex')
                title('Obj. function evolution','Interpreter','latex')
                set(gca,'TickLabelInterpreter','latex')
                box on
                subplot(1,3,3)
                hold on
                plot(obj.iterVec,obj.cnstVec1,'b',obj.iterVec,obj.cnstVec2,'g')
                hold off
                xlabel('Iteration','Interpreter','latex')
                ylabel('Constraint violation','Interpreter','latex')
                title('Constr. violation evolution','Interpreter','latex')
                legend('Constraint 1', 'Constraint 2','Interpreter','latex')
                set(gca,'TickLabelInterpreter','latex')
                box on
                drawnow
            else
                obj.xVec1 = [obj.xVec1;x(1)];
                obj.iterVec = [obj.iterVec;iter];
                obj.cVec = [obj.cVec;c];
                obj.cnstVec1 = [obj.cnstVec1;cnst(1)];
                subplot(1,3,1)
                hold on
                plot(obj.xVec1,'r')
                hold off
                xlabel('$x_1$','Interpreter','latex')
                title('Design variables','Interpreter','latex')
                set(gca,'TickLabelInterpreter','latex')
                subplot(1,3,2)
                hold on
                plot(obj.iterVec,obj.cVec,'b')
                hold off
                xlabel('Iteration','Interpreter','latex')
                ylabel('Objective function $J(x)$','Interpreter','latex')
                title('Obj. function evolution','Interpreter','latex')
                set(gca,'TickLabelInterpreter','latex')
                subplot(1,3,3)
                hold on
                plot(obj.iterVec,obj.cnstVec1,'b')
                hold off
                xlabel('Iteration','Interpreter','latex')
                ylabel('Constraint violation','Interpreter','latex')
                title('Constr. violation evolution','Interpreter','latex')
                legend('Constraint','Interpreter','latex')
                set(gca,'TickLabelInterpreter','latex')
                drawnow
            end
            vector1 = obj.xVec1;
            vector2 = obj.xVec2;
            iterations = obj.iterVec;
            costVec = obj.cVec;
            constrVec1 = obj.cnstVec1;
            constrVec2 = obj.cnstVec2;
%             save('Cost_Const_Test4.mat','iterations',"costVec","constrVec1","constrVec2")
%             save('DesignVariableEvolution_Test1.mat','vector1','vector2','-mat')
%             obj.xVec1 = [obj.xVec1;x(1)];
%             obj.xVec2 = [obj.xVec2;x(2)];
%             if iter == 0
%                 open FeasibleRegionTest1_2.fig
%             end
%             figure(10)
%             hold on
%             v = 0:0.001:5;
%             [x y] = meshgrid(v);
%             cond1 = 1./x - y < 0;
%             cond2 = x+ y - 3 < 0;
%             cond1 = double(cond1);
%             cond2 = double(cond2);
%             cond1(cond1 == 0) = NaN;
%             cond2(cond2 == 0) = NaN;
%             cond = cond1.*cond2;
%             h = surf(x,y,cond);
%             colormap([0,1,1])
%             set(h,'edgecolor','none')
%             view(0,90)
%             plot(obj.xVec1,obj.xVec2)
%             hold off

        end
                
    end
    
end