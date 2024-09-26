classdef LagrangianPlotter < handle

    properties (Access = private)
        lagrangianFunction
    end

    methods (Access = public)

        function obj = LagrangianPlotter(cParams)
            obj.init(cParams)
        end

        function plot(obj)
            lf = obj.lagrangianFunction;
            coord  = lf.getCoord();
            connec = lf.getConnec();
            zv     = lf.fValues;
            ndimf  = lf.ndimf;
            
            switch lf.getOrderInText()
                case 'CONSTANT'
                    f = lf.project('P1D');
                    f.plot()

                case 'LINEAR'
                    switch lf.mesh.type
                        case {'TRIANGLE','QUAD'}
                            x = coord(:,1);
                            y = coord(:,2);
                            figure()
                            for idim = 1:ndimf
                                subplot(1,ndimf,idim);
                                z = zv(:,idim);
                                a = trisurf(connec,x,y,z);
                                view(0,90)
                                %             colorbar
                                shading interp
                                a.EdgeColor = [0 0 0];
                                title(['dim = ', num2str(idim)]);
                            end
                        case 'LINE'
                            coord = coord(:,1);
                            z = zv;
                            figure()
                            plot(coord,z)
                    end

                otherwise
                    figure()
                    for idim = 1:ndimf
                        subplot(1,ndimf,idim);
                        hold on
                        x = coord(:,1);
                        y = coord(:,2);
                        z = zv(:,idim);
                        %better to remesh (now only plotting the linear part)
                        T = lf.mesh.connec;
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
            obj.lagrangianFunction = cParams.function;
        end

    end

end