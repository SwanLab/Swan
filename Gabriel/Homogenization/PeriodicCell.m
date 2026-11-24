classdef PeriodicCell < handle

    properties (GetAccess = public, SetAccess = private)
        mesh
        masterSlave
        vertexNodes
        divisions
        basis
    end

    properties (Constant, Access = private)
        vertexTolerance    = 1e-10
        reductionTolerance = 1e-12
    end

    methods (Access = public)

        function obj = PeriodicCell(cParams)
            [obj.basis,obj.divisions] = PeriodicCell.discretize(cParams.T,cParams.divUnit,cParams.strategy);
            obj.createMesh();
        end

    end

    methods (Access = public, Static)

        function [B,div] = discretize(T,divUnit,strategy)
            switch strategy
                case 'Mapped'
                    B   = T;
                    div = PeriodicCell.toEven(divUnit*[1,1]);
                case 'Adapted'
                    B   = T;
                    div = PeriodicCell.toEven(divUnit*PeriodicCell.computeAdaptedRatio(T));
                case 'Reduced'
                    B   = PeriodicCell.reduceBasis(T);
                    div = PeriodicCell.toEven(divUnit*vecnorm(B)/sqrt(abs(det(B))));
                otherwise
                    error('PeriodicCell:UnknownStrategy','Unknown mesh strategy: %s',strategy);
            end
        end

    end

    methods (Access = private)

        function createMesh(obj)
            s.T        = obj.basis;
            s.div      = obj.divisions;
            s.filename = '';
            mc = MeshCreator(s);
            mc.computeMeshNodes();
            sM.coord        = mc.coord;
            sM.connec       = mc.connec;
            obj.mesh        = Mesh.create(sM);
            obj.masterSlave = mc.masterSlaveIndex;
            obj.vertexNodes = obj.findVertexNodes((obj.basis\mc.coord')');
        end

        function nodes = findVertexNodes(obj,xi)
            distanceToCorner = min(abs(xi),abs(xi-1));
            nodes = find(all(distanceToCorner < obj.vertexTolerance,2));
        end

    end

    methods (Access = private, Static)

        function B = reduceBasis(T)
            B = T;
            while true
                if norm(B(:,1)) > norm(B(:,2))
                    B = B(:,[2 1]);
                end
                mu = dot(B(:,1),B(:,2))/dot(B(:,1),B(:,1));
                if abs(mu) <= 0.5 + PeriodicCell.reductionTolerance
                    break
                end
                B(:,2) = B(:,2) - round(mu)*B(:,1);
            end
            if det(B) < 0
                B(:,2) = -B(:,2);
            end
        end

        function ratio = computeAdaptedRatio(T)
            lengths  = vecnorm(T);
            cosAngle = abs(dot(T(:,1),T(:,2)))/prod(lengths);
            ratio    = max(1,(lengths/min(lengths))*cosAngle);
        end

        function n = toEven(x)
            n = 2*ceil(x/2 - 1e-9);
        end

    end

end