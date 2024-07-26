classdef Display_Surf < Display_Abstract

    properties (Access = private)
        mesh
        barLim
    end

    methods (Access = public)

        function obj = Display_Surf(cParams)
            obj@Display_Abstract(cParams.title)
            obj.mesh = cParams.mesh;
            if isfield(cParams,'barLim')
                obj.barLim = cParams.barLim;
            else
                obj.barLim = [];
            end
        end
    end

    methods (Access = protected)
        
        function setChartType(obj)
            z = zeros(obj.mesh.nnodes,1);
            a = obj.createTrisurf(z);
            obj.handle = a;
        end
        
    end

    methods (Access = public)

        function updateParams(obj,it,value)
            obj.iterationArray = it;
            if ~isempty(value)
                obj.valueArray = value;
            end
        end

        function refresh(obj)
            if ~isempty(obj.valueArray) && ~isempty(obj.iterationArray)
                z = obj.valueArray;
                a = obj.createTrisurf(z);
                obj.handle = a;
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

            if ~isempty(obj.iterationArray)
                iter = num2str(obj.iterationArray);
                title(strcat(obj.figTitle,' iter:',iter));
            else
                title(obj.figTitle);
            end
        end

    end

end