classdef DisplaySurf < DisplayAbstract

    properties (Access = private)
        fun
        barLim
        faces

        FieldData
        iter
    end

    methods (Access = public)

        function obj = DisplaySurf(cParams)
            obj@DisplayAbstract(cParams)
            obj.fun = cParams.fun;
            if ~isempty(cParams.barLim)
                obj.barLim = cParams.barLim;
            else
                obj.barLim = "auto";
            end
        end
    end

    methods (Access = protected)
        
        function setChartType(obj)
            obj.FieldData = zeros(size(obj.fun.fValues(:,1)));
            a = obj.createTrisurf(obj.FieldData);
            obj.handle = a;
            obj.faces = obj.handle.Faces;
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
                set(obj.handle,'ZData',obj.FieldData,'CData',obj.FieldData,'Faces',obj.faces);
            end
        end
    end

    methods (Access = private)

        function a = createTrisurf(obj,z)
            connec = obj.fun.mesh.connec;
            x = obj.fun.mesh.coord(:,1);
            y = obj.fun.mesh.coord(:,2);
            
            a = trisurf(connec,x,y,z);
            view(0,90)
            shading interp
            a.EdgeColor = [0 0 0];
            colorbar
            clim(obj.barLim);
            hold off
        end

    end

end