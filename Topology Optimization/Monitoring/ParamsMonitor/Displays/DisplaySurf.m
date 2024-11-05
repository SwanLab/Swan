classdef DisplaySurf < DisplayAbstract

    properties (Access = private)
        mesh
        barLim
        faces

        FieldData
        iter
    end

    methods (Access = public)

        function obj = DisplaySurf(cParams)
            obj@DisplayAbstract(cParams.title,cParams.position)
            obj.mesh = cParams.mesh;
            if isfield(cParams,"barLim")
                obj.barLim = cParams.barLim;
            else
                obj.barLim = "auto";
            end
        end
    end

    methods (Access = protected)
        
        function setChartType(obj)
            obj.FieldData = zeros(obj.mesh.nnodes,1);
            a = obj.createTrisurf(obj.FieldData);
            obj.handle = a;
            obj.faces = obj.handle.Faces;
        end
        
    end

    methods (Access = public)

        function updateParams(obj,it,value) 
            if ~isempty(value)
                obj.FieldData = value;
                obj.iter = it;
            end
        end

        function refresh(obj)
            if ~isempty(obj.FieldData) && ~isempty(obj.iter)
                t = strcat(obj.figTitle,' iter:',num2str(obj.iter));
                axis = findobj(gcf,'Type','Axes');
                set(axis(obj.position).Title,'String',t);
                set(obj.handle,'ZData',obj.FieldData,'CData',obj.FieldData,'Faces',obj.faces);
                drawnow
            end
        end
    end

    methods (Access = private)

        function a = createTrisurf(obj,z)
            connec = obj.mesh.connec;
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            
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