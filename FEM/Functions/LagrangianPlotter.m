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
            zv     = obj.fValues;
            
            switch lf.getOrderInText()
                case 'CONSTANT'
                    f = lf.project('P1D');
                    f.plot()

                case 'LINEAR'
                    switch lf.mesh.type
                        case {'TRIANGLE','QUAD'}
                            x = obj.coord(:,1);
                            y = obj.coord(:,2);
                            figure()
                            for idim = 1:obj.ndimf
                                subplot(1,obj.ndimf,idim);
                                z = zv(:,idim);
                                a = trisurf(obj.connec,x,y,z);
                                view(0,90)
                                %             colorbar
                                shading interp
                                a.EdgeColor = [0 0 0];
                                title(['dim = ', num2str(idim)]);
                            end
                        case 'LINE'
                            x = obj.coord(:,1);
                            z = zv;
                            figure()
                            plot(x,z)
                    end

                otherwise
                    figure()
                    for idim = 1:obj.ndimf
                        subplot(1,obj.ndimf,idim);
                        hold on
                        x = obj.coord(:,1);
                        y = obj.coord(:,2);
                        z = zv(:,idim);
                        %better to remesh (now only plotting the linear part)
                        T = obj.connec;
                        a = trisurf(T,x,y,z);
                        view(0,90)
                        colorbar
                        shading interp
                        grid on
                        title(['dim = ', num2str(idim)]);
                        a.EdgeColor = [0 0 0];
                        a.EdgeColor = [0 0 0];
                    end
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