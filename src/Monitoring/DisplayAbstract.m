classdef DisplayAbstract < handle
    
    properties (Access = protected)
        figTitle

        position
        figureIdx
        handle
        iteration
        ArrayDataX
        ArrayDataY
        style
    end

    methods (Access = protected, Abstract)
        setChartType(obj)
    end

    methods (Access = public, Abstract)
        updateParams(ibj,it,value)
    end

    methods (Access = public , Static)
        
        function obj = create(s)
            d = DisplayFactory();
            obj = d.create(s);
        end
    end

    methods (Access = public)

        function obj = DisplayAbstract(cParams)
            obj.figTitle = cParams.title;
            obj.position = cParams.position;
            obj.figureIdx = cParams.figureIdx;
        end

        function show(obj,nRows,nCols,i,margins)
            obj.setStyle(nRows,nCols,i,margins);
            obj.setChartType();
            title(obj.figTitle);
            grid on
        end

    end

    methods (Access = public, Abstract)
        refresh(obj)
    end

    methods (Access = protected)

        function DispAxes = obtainDisplayAxes(obj)
            AllFigs = findobj('Type','Figure');
            Fig = findobj(AllFigs,'Number',obj.figureIdx);
            axes = flip(findobj(Fig,'Type','axes'));
            DispAxes = axes(obj.position);
        end

    end
    
    methods (Access = private)
        
        function setStyle(obj,nRows,nCols,i,margins)
            obj.style = subplot_tight(nRows,nCols,i,margins);
        end

    end
    
end


