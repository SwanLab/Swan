classdef AcademicProblemPathComputer < handle
    
    properties (Access = private)
        printingPath
        cost
        constraint
        constraintCase
        boxConstraints
    end

    properties (Access = private)
        colors = 'mgr'
    end

    methods (Access = public)
        function obj = AcademicProblemPathComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj,values)
            if obj.printingPath && size(values,1) == 2 % Just f(x,y) problems
                x     = values(1,:);
                y     = values(2,:);
                [X,Y] = obj.setMeshGrid(x,y);
                obj.plotCostContour(X,Y);
                isF = obj.plotConstraintContour(X,Y);
                obj.plotBoxConstraints(X,Y,isF);
                obj.plotDesignVariablePath(values);
            end
        end
    end

    methods (Access = private)
        function init(obj, cParams)
            obj.printingPath      = cParams.printingPath;
            obj.cost              = cParams.cost.cF;
            obj.constraint        = cParams.constraint.cF;
            obj.constraintCase    = cParams.settings.constraintCase;
            obj.boxConstraints.ub = cParams.ub;
            obj.boxConstraints.lb = cParams.lb;
        end

        function plotCostContour(obj,X,Y)
            x  = @(i) X.*(i==1)+Y.*(i==2);
            fX = obj.cost(x);
            v  = linspace(min(fX(:)),max(fX(:)),30);
            figure
            contour(X,Y,fX,v,'linewidth',1);
            hold on;
        end

        function isFeasible = plotConstraintContour(obj,X,Y)
            x1 = unique(X);
            x2 = unique(Y);
            x  = @(i) X.*(i==1)+Y.*(i==2);
            nConst = length(obj.constraint);
            isFeasible = ones(length(x1),length(x2));
            for i = 1:nConst
                c  = obj.constraint{i};
                cr = obj.colors(i);
                fX = c(x);
                v  = linspace(-1e-3,1e-3,2);
                [~,h] = contour(X,Y,fX,v,'linewidth',2);
                h.LineColor=cr;
                hold on;
                switch obj.constraintCase{i}
                    case 'INEQUALITY'
                        isFeasible = isFeasible & fX<=0;
                end
            end
        end

        function plotBoxConstraints(obj,X,Y,isFeasible)
            x1 = unique(X);
            x2 = unique(Y);
            x  = @(i) X.*(i==1)+Y.*(i==2);
            v  = linspace(-1e-3,1e-3,2);

            xlb = obj.boxConstraints.lb(1);
            xub = obj.boxConstraints.ub(1);
            ylb = obj.boxConstraints.lb(2);
            yub = obj.boxConstraints.ub(2);

            c{1}  = @(x) xlb-x(1);
            c{2}  = @(x) x(1)-xub;
            c{3}  = @(x) ylb-x(2);
            c{4}  = @(x) x(2)-yub;

            for i = 1:length(c)
                fX = c{i}(x);
                [~,h] = contour(X,Y,fX,v,'linewidth',1);
                h.LineColor='r';
                hold on;
                isFeasible = isFeasible & fX<=0;
            end

            output = NaN(length(x1),length(x2));
            output(isFeasible) = 1;
            a = surf(X,Y,output,'EdgeColor','none');
            alpha(.25);
            view([0 90]);
            set(get(get(a,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            hold on;
        end

        function plotDesignVariablePath(obj,values)
            vx = values(1,:);
            vy = values(2,:);
            plot(vx(1),vy(1),"o",'MarkerSize',10,'MarkerFaceColor','red');
            iters = length(vx);
            dit   = floor(iters/11);
            for i = 1:10
                plot(vx(i*dit),vy(i*dit),"o",'MarkerSize',5,'MarkerFaceColor','red');
            end
            plot(vx,vy,'-k','linewidth',1);
            plot(vx(end),vy(end),"p",'MarkerSize',10,'MarkerFaceColor','red');
            nConst = length(obj.constraint);
            leg = "J";
            for i = 1:nConst
                leg = [leg,string(['g_',char(string(i))])];
            end
            leg = [leg,"Box","","",""];
            leg = [leg,"Initial point","Intermediate iter","","","","","","","","","","Path","Solution"];
            legend(leg);
            hold off;
        end
    end

    methods (Static, Access = private)
        function [X,Y] = setMeshGrid(x,y)
            xmin  = min(x);
            xmax  = max(x);
            dx    = xmax-xmin;
            ymin  = min(y);
            ymax  = max(y);
            dy    = ymax-ymin;
            xV    = linspace(-3.05,6.05,2000);
            yV    = linspace(-1.05,3.05,2000);
            [X,Y] = meshgrid(xV,yV);
        end
    end
end