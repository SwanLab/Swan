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
                obj.plotConstraintContour(x);
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
            contour(X,Y,fX,v,'linewidth',2);
            hold on;
        end

        function plotConstraintContour(obj,x)
            xmin   = min(x);
            xmax   = max(x);
            xV     = linspace(0.8*xmin,1.2*xmax,2000);
            x      = @(y,i) xV.*(i==1)+y.*(i==2);
            nConst = length(obj.constraint);
            for i = 1:nConst
                c   = obj.constraint{i};
                fun = @(y) c(x(y,':'));
                y   = fzero(fun,0);
                % here
            end
        end
    end

    methods (Static, Access = private)
        function [X,Y] = setMeshGrid(x,y)
            xmin  = min(x);
            xmax  = max(x);
            ymin  = min(y);
            ymax  = max(y);
            xV    = linspace(0.8*xmin,1.2*xmax,2000);
            yV    = linspace(0.8*ymin,1.2*ymax,2000);
            [X,Y] = meshgrid(xV,yV);
        end
    end

end