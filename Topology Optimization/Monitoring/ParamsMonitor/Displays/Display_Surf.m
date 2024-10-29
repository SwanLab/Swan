classdef Display_Surf < Display_Abstract

    properties (Access = private)
        mesh
        faces
        vertices
        barLim
    end

    methods (Access = public)

        function obj = Display_Surf(cParams)
            obj@Display_Abstract(cParams.title)
            obj.mesh = cParams.mesh;
            obj.barLim = cParams.barLim;
        end
    end

    methods (Access = protected)
        
        function setChartType(obj)
            z = zeros(obj.mesh.nnodes,1);
            a = obj.createTrisurf(z);
            obj.handle = a;
            obj.faces = obj.handle.Faces;
            obj.vertices = obj.handle.Vertices;
        end
        
    end

    methods (Access = public)

        function updateParams(obj,it,value) 
            if ~isempty(value)
                obj.valueArray = value;
                obj.iterationArray = it;
            end
        end

        function refresh(obj)
            if ~isempty(obj.valueArray) && ~isempty(obj.iterationArray)
                z = obj.valueArray;
                obj.vertices(:,3) = z;
                set(obj.handle,'ZData',z,'CData',z,'Faces',obj.faces,'Vertices',obj.vertices);
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

            if ~isempty(obj.iterationArray)
                iter = num2str(obj.iterationArray);
                title(strcat(obj.figTitle,' iter:',iter));
            else
                title(obj.figTitle);
            end
        end

    end

end