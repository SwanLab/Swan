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
            figure()
            for idim = 1:obj.ndimf
                subplot(1,obj.ndimf,idim);
                x  = obj.coord{idim}(:,1);
                y  = obj.coord{idim}(:,2);  
                z  = obj.fValues(:,idim);            
                a = trisurf(obj.connec{idim},x,y,z);
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