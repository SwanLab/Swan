classdef AcademicMonitoring < handle

    properties (Access = private)
        designVarVec
        iterVec
        costVec
        dualVariableVec
        cnstVec
    end

    properties (Access = private)
        cost
        constraint
        designVariable
        dualVariable
        shallPrint
        colorsSample
    end
    
    methods (Access = public)
        function obj = AcademicMonitoring(cParams)
            obj.init(cParams);
        end

        function create(obj,cParams)
            
        end

        function compute(obj,cParams)
            if (obj.shallPrint)
                obj.plot(cParams);
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.cost           = cParams.cost;
            obj.constraint     = cParams.constraint;
            obj.designVariable = cParams.designVar;
            obj.dualVariable   = cParams.dualVariable;
            obj.shallPrint     = cParams.shallPrint;
            obj.colorsSample   = ["r","b","g","m","k","c","y"];
        end
        
        function plot(obj,cParams)
            obj.updateNewIterationValues(cParams)
            obj.plotDesignVariables();
            obj.plotCostFunction();
            obj.plotDualVariable();
            obj.plotConstraints();
            drawnow
        end

        function updateNewIterationValues(obj,cParams)
            iter                = cParams.nIter;
            x                   = obj.designVariable.value;
            c                   = obj.cost.value;
            l                   = obj.dualVariable.value;
            cnst                = obj.constraint.value;
            obj.designVarVec    = [obj.designVarVec,x];
            obj.iterVec         = [obj.iterVec,iter];
            obj.costVec         = [obj.costVec,c];
            obj.dualVariableVec = [obj.dualVariableVec,l'];
            obj.cnstVec         = [obj.cnstVec,cnst];
        end

        function plotDesignVariables(obj)
            subplot(1,4,1)
            hold on
            dim           = size(obj.designVariable.value,1);
            legendContent = strings(1,dim);
            for i = 1:dim
                plot(obj.iterVec,obj.designVarVec(i,:),obj.colorsSample(i))
                legendContent(i) = string(['x',num2str(i)]);
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
            plot(obj.iterVec,obj.costVec,'b')
            hold off
            xlabel('Iteration')
            ylabel('Objective function J(x)')
            title('Objective function evolution')
        end

        function plotDualVariable(obj)
            subplot(1,4,3)
            hold on
            dim           = obj.constraint.nSF;
            legendContent = strings(1,dim);
            for i = 1:dim
                plot(obj.iterVec,obj.dualVariableVec(i,:),obj.colorsSample(i))
                legendContent(i) = string(['$\lambda$ ',num2str(i)]);
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
            dim           = obj.constraint.nSF;
            legendContent = strings(1,dim);
            for i = 1:dim
                plot(obj.iterVec,obj.cnstVec(i,:),obj.colorsSample(i))
                legendContent(i) = string(['Constraint ',num2str(i)]);
            end
            hold off
            xlabel('Iteration')
            ylabel('Constraint violation')
            title('Constraints violation evolution')
            legend(legendContent)
        end    
    end
end