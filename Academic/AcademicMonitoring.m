classdef AcademicMonitoring < handle

    properties (Access = private)
        dimDesignVar
        dimConstr
        colorsSample
        xVec
        iterVec
        cVec
        dualVariableVec
        cnstVec
    end

    properties (Access = private)
        cost
        constraint
        designVariable
        dualVariable
        shallPrint
    end
    
    methods (Access = public)
        
        function obj = AcademicMonitoring(cParams)
            obj.init(cParams);
        end

        function create(obj,cParams)
            x                = cParams.designVar.value;
            obj.dimDesignVar = size(x,1);
            obj.dimConstr    = cParams.nConstraints;
            obj.colorsSample = ["r","b","g","m","k","c","y"];
        end

        function compute(obj,cParams)
            if (obj.shallPrint)
                obj.plot(cParams);
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.cost              = cParams.cost;
            obj.constraint        = cParams.constraint;
            obj.designVariable    = cParams.designVar;
            obj.dualVariable      = cParams.dualVariable;
            obj.shallPrint        = cParams.shallPrint;
        end
        
        function plot(obj,cParams)
            iter        = cParams.nIter;
            x           = obj.designVariable.value;
            c           = obj.cost.value;
            l           = obj.dualVariable.value;
            cnst        = obj.constraint.value;
            obj.xVec    = [obj.xVec,x];
            obj.iterVec = [obj.iterVec,iter];
            obj.cVec    = [obj.cVec,c];
            obj.dualVariableVec = [obj.dualVariableVec,l'];
            obj.cnstVec = [obj.cnstVec,cnst];
            obj.plotDesignVariables();
            obj.plotCostFunction();
            obj.plotDualVariable();
            obj.plotConstraints();
            drawnow
        end

        function plotDesignVariables(obj)
            subplot(1,4,1)
            hold on
            legendContent = [];
            for i = 1:obj.dimDesignVar
                plot(obj.iterVec,obj.xVec(i,:),obj.colorsSample(i))
                legendContent = [legendContent,string(['x',num2str(i)])];
            end
            hold off
            xlabel('Iteration')
            ylabel('Design variable')
            title('Design variables evolution')
            legend(legendContent)
        end

        function plotCostFunction(obj)
            subplot(1,4,2)
            hold on
            plot(obj.iterVec,obj.cVec,'b')
            hold off
            xlabel('Iteration')
            ylabel('Objective function J(x)')
            title('Objective function evolution')
        end

        function plotDualVariable(obj)
            subplot(1,4,3)
            hold on
            legendContent = [];
            for i = 1:obj.dimConstr
                plot(obj.iterVec,obj.dualVariableVec(i,:),obj.colorsSample(i))
                legendContent = [legendContent,string(['$\lambda$ ',num2str(i)])];
            end
            hold off
            xlabel('Iteration')
            ylabel('Constraint violation')
            title('Constraints violation evolution')
            legend(legendContent,'interpreter','latex')
        end

        function plotConstraints(obj)
            subplot(1,4,4)
            hold on
            legendContent = [];
            for i = 1:obj.dimConstr
                plot(obj.iterVec,obj.cnstVec(i,:),obj.colorsSample(i))
                legendContent = [legendContent,string(['Constraint ',num2str(i)])];
            end
            hold off
            xlabel('Iteration')
            ylabel('Constraint violation')
            title('Constraints violation evolution')
            legend(legendContent)
        end
                
    end
    
end