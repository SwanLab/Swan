classdef DisplayQuiver < DisplayAbstract

    properties (Access = private)
        fun
        FieldData
        iter
    end

    methods (Access = public)

        function obj = DisplayQuiver(cParams)
            obj@DisplayAbstract(cParams)
            obj.fun = cParams.fun;
        end
    end

    methods (Access = protected)
        
        function setChartType(obj)
            nnodes = obj.fun.nDofs/obj.fun.ndimf;
            obj.FieldData = zeros(nnodes,2);
            a = obj.createQuiver(obj.FieldData);
            obj.handle = a;
        end

    end

    methods (Access = public)

        function updateParams(obj,it,fValues) 
                obj.FieldData = fValues;
                obj.iter = it;
        end

        function refresh(obj)
            if ~isempty(obj.FieldData) && ~isempty(obj.iter)
                t = strcat(obj.figTitle,' / iter:',num2str(obj.iter));
                axes = obj.obtainDisplayAxes();
                set(axes.Title,'String',t);
                set(obj.handle,'UData',obj.FieldData(:,1),'VData',obj.FieldData(:,2));
                shading(axes,"interp")
            end
        end
    end

    methods (Access = private)

        function a = createQuiver(obj,z)
            n = 1;
            coord = obj.fun.getDofCoord; 
            x = coord(1:n:end,1);
            y = coord(1:n:end,2);
            fX = z(1:n:end,1);
            fY = z(1:n:end,2);
            a = quiver(x, y, fX, fY, 'AutoScale', 'on', 'LineWidth', 1.5);              
            axis equal;  
            box on;     
            xlim([min(x)-1, max(x)+1]);
            ylim([min(y)-1, max(y)+1]); 
        end

    end

end