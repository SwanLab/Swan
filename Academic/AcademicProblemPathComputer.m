classdef AcademicProblemPathComputer < handle
    
    properties (Access = private)
        printingPath
        cost
        constraint
        constraintCase
    end

    methods (Access = public)
        function obj = AcademicProblemPathComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj,values)
            if obj.printingPath && size(values,1) == 2
                x     = values(1,:);
                y     = values(2,:);
                [X,Y] = obj.setMeshGrid(x,y);
                obj.plotCostContour(X,Y);
                obj.plotConstraintContour(X,Y);
                obj.plotDesignVariablePath(values);
            end
        end
    end

    methods (Access = private)
        function init(obj, cParams)
            obj.printingPath   = cParams.printingPath;
            obj.cost           = cParams.cost.cF;
            obj.constraint     = cParams.constraint.cF;
            obj.constraintCase = cParams.settings.constraintCase;
        end

        function plotCostContour(obj,X,Y)
            x  = @(i) X.*(i==1)+Y.*(i==2);
            fX = obj.cost(x);
            v  = linspace(min(fX(:)),max(fX(:)),30);
            figure
            contour(X,Y,fX,v,'linewidth',1);
            hold on;
        end

        function plotConstraintContour(obj,X,Y)
            x  = @(i) X.*(i==1)+Y.*(i==2);
            nConst = length(obj.constraint);
            for i = 1:nConst
                c  = obj.constraint{i};
                fX = c(x);
                v  = linspace(-1e-3,1e-3,2);
                [C,h] = contour(X,Y,fX,v,'linewidth',2);
                h.LineColor='m';
                hold on;
                switch obj.constraintCase{i}
                    case 'INEQUALITY'
                        ymin = min(Y(:));
                        n  = size(C,2);
                        xg = C(1,2:n/2);
                        yg = C(2,2:n/2);
                        area(xg, yg,BaseValue=ymin,FaceColor='m');
                        alpha(.4);
                end
            end
        end

        function plotDesignVariablePath(obj,values)
            vx = values(1,:);
            vy = values(2,:);
            plot(vx(1),vy(1),"o",'MarkerSize',10,'MarkerFaceColor','red');
            plot(vx(2:end-1),vy(2:end-1),'-k','linewidth',1);
            plot(vx(end),vy(end),"p",'MarkerSize',10,'MarkerFaceColor','red');
            nConst = length(obj.constraint);
            leg = "J";
            for i = 1:nConst
                leg = [leg,string(['g_',char(string(i))])];
            end
            leg = [leg,"Initial point","Path","Solution"];
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
            xV    = linspace(xmin-0.1*dx,xmax+0.1*dx,2000);
            yV    = linspace(ymin-0.1*dy,ymax+0.1*dy,2000);
            [X,Y] = meshgrid(xV,yV);
        end
    end
end