classdef ConsistentLatticeVolumeConstraint < handle
    properties (Access = private)
        mesh
        volumeTarget
        sideFields
        totalVolume
        consistentMass
        lumpedMass
    end
    methods (Access = public)
        function obj = ConsistentLatticeVolumeConstraint(cParams)
            obj.init(cParams);
        end
        function [J,dJ] = computeFunctionAndGradient(obj,x)
            xD = x.obtainDomainFunction();
            m1 = xD{obj.sideFields(1)};
            m2 = xD{obj.sideFields(2)};
            V  = obj.computeVolumeFraction(m1,m2);
            J  = V/obj.volumeTarget - 1;
            dJ = obj.computeGradient(xD,m1,m2);
        end
    end
    methods (Access = private)
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.volumeTarget = cParams.volumeTarget;
            obj.sideFields   = cParams.sideFields;
            obj.totalVolume  = obj.mesh.computeVolume();
            obj.createMassMatrices();
        end
        function V = computeVolumeFraction(obj,m1,m2)
            hole = m1.*m2;
            V    = 1 - Integrator.compute(hole,obj.mesh,2)/obj.totalVolume;
        end
        function dJ = computeGradient(obj,xD,m1,m2)
            scale = obj.volumeTarget*obj.totalVolume;
            dJ    = cell(numel(xD),1);
            for i = 1:numel(xD)
                dJ{i} = copy(xD{i});
                dJ{i}.setFValues(zeros(size(xD{i}.fValues)));
            end
            dJ{obj.sideFields(1)}.setFValues(-obj.project(m2.fValues)/scale);
            dJ{obj.sideFields(2)}.setFValues(-obj.project(m1.fValues)/scale);
        end
        function p = project(obj,values)
            p = (obj.consistentMass*values(:))./obj.lumpedMass;
        end
        function createMassMatrices(obj)
            coord  = obj.mesh.coord;
            connec = obj.mesh.connec;
            nNodes = size(coord,1);
            e1     = coord(connec(:,2),:) - coord(connec(:,1),:);
            e2     = coord(connec(:,3),:) - coord(connec(:,1),:);
            area   = abs(e1(:,1).*e2(:,2) - e1(:,2).*e2(:,1))/2;
            M      = sparse(nNodes,nNodes);
            for i = 1:3
                for j = 1:3
                    M = M + sparse(connec(:,i),connec(:,j),area*(1+(i==j))/12,nNodes,nNodes);
                end
            end
            obj.consistentMass = M;
            obj.lumpedMass     = full(sum(M,2));
        end
    end
    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume constraint';
        end
    end
end