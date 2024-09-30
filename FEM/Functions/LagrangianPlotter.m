classdef LagrangianPlotter < handle

    properties (Access = private)
        coord
        connec
        fValues
        ndimf
    end

    methods (Access = public)

        function obj = LagrangianPlotter(cParams)
            obj.init(cParams)
        end

        function plot(obj)
            x  = obj.coord(:,1);
            y  = obj.coord(:,2);
            z  = obj.fValues;            
            figure()
            for idim = 1:obj.ndimf
                subplot(1,obj.ndimf,idim);
                zi = z(:,idim);
                a = trisurf(obj.connec,x,y,zi);
                view(0,90)
                %colorbar
                shading interp
                a.EdgeColor = [0 0 0];
                title(['dim = ', num2str(idim)]);
            end       
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.coord   = cParams.coord;
            obj.connec  = cParams.connec;
            obj.fValues = cParams.fValues;
            obj.ndimf   = cParams.ndimf;
        end

    end

end