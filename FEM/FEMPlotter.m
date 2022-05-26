classdef FEMPlotter < handle

    properties (Access = private)
        dim
        mesh
        displacement
        dispCoords
    end

    methods (Access = public)

        function obj = FEMPlotter(cParams)
            obj.init(cParams);
        end

        function plot(obj)
            ndimf = obj.dim.ndimField;
            switch ndimf
                case 2
                    obj.plotFem2D();
                case 3
                    obj.plotFem3D();
            end
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.dim          = cParams.dim;
            obj.mesh         = cParams.mesh;
            obj.displacement = cParams.displacement;
            obj.dispCoords   = obj.computeDisplacedCoords();
        end

        function dispCoords = computeDisplacedCoords(obj)
            ndimf = obj.dim.ndimField;
            ndof = obj.dim.ndof;
            coords = zeros(ndof,1);
            for i = 1:ndimf
                dofs = i:ndimf:ndof;
                coor = obj.mesh.coord(:,i);
                coords(dofs) = coor;
            end
            delta = obj.displacement;
            dispCoords = coords + delta;
        end

        function plotFem2D(obj)
            obj.plotNodes2D();
            obj.plotDisplacement2D();
        end

        function plotFem3D(obj)
            obj.plotNodes3D();
            obj.plotDisplacement3D();
        end

        function plotNodes2D(obj)
            Tn = obj.mesh.connec;
            x  = obj.mesh.coord(:,1);
            y  = obj.mesh.coord(:,2);
            figure()
            hold on
            colormap jet;
            plot(x(Tn)',y(Tn)','--','linewidth',0.5);
        end

        function plotDisplacement2D(obj)
            ndimf  = obj.dim.ndimField;
            ndof   = obj.dim.ndof;
            coords = obj.dispCoords;
            x = coords(1:ndimf:ndof);
            y = coords(2:ndimf:ndof);
            Tn   = obj.mesh.connec;
            plot(x(Tn)',y(Tn)','-k','linewidth',0.5);
        end

        function plotNodes3D(obj)
            Tn = obj.mesh.connec;
            x  = obj.mesh.coord(:,1);
            y  = obj.mesh.coord(:,2);
            z  = obj.mesh.coord(:,3);
            figure()
            hold on
            colormap jet;
            plot3(x(Tn)',y(Tn)',z(Tn)','--','linewidth',0.5);
        end

        function plotDisplacement3D(obj)
            ndimf  = obj.dim.ndimf;
            ndof   = obj.dim.ndof;
            coords = obj.dispCoords;
            x = coords(1:ndimf:ndof);
            y = coords(2:ndimf:ndof);
            z = coords(3:ndimf:ndof);
            Tn   = obj.mesh.connec;
            plot3(x(Tn)',y(Tn)',z(Tn)','-k','linewidth',0.5);
        end

    end
end