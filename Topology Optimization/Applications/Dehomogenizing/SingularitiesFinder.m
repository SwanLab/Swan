classdef SingularitiesFinder < handle

    properties (Access = public)

    end

    properties (Access = private)
        orientation
        mesh
    end

    properties (Access = private)
        meshDisc
        isElemSingular
    end

    methods (Access = public)

        function obj = SingularitiesFinder(cParams)
            obj.init(cParams)
        end

        function isS = computeSingularElements(obj)
            obj.computeSingularities();
            isS = obj.isElemSingular;
        end

        function plot(obj)
            obj.createDiscontinousMesh();
            obj.plotOrientationVector();
            obj.plotSingularities();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh        = cParams.mesh;
            obj.orientation = cParams.orientation;
        end

        function createDiscontinousMesh(obj)
            m = obj.mesh.createDiscontinousMesh();
            obj.meshDisc = m;
        end

        function computeSingularities(obj)
            b = obj.computeOrientationOfDoubleAngle;
            b = obj.mapP1ToP1Discontinous(b);
            b1 = b(:,:,1);
            b2 = b(:,:,2);
            b3 = b(:,:,3);
            b1b2 = obj.scalarProduct(b1,b2);
            b1b3 = obj.scalarProduct(b1,b3);
            b2b3 = obj.scalarProduct(b2,b3);
            isS = sign(b1b2.*b1b3.*b2b3);
            obj.isElemSingular = isS;
        end

        function b = computeOrientationOfDoubleAngle(obj)
            beta = obj.computeBetaAngle();
            b(:,1) = cos(beta);
            b(:,2) = sin(beta);
        end

        function beta = computeBetaAngle(obj)
            a = obj.orientation;            
            alpha = atan2(a(:,2),a(:,1));
            beta = 2*alpha;
        end

        function plotOrientationVector(obj)
            figure()
            a = obj.orientation;
            x = obj.mesh.coord(:,1);
            y = obj.mesh.coord(:,2);
            ax = a(:,1);
            ay = a(:,2);
            quiver(x,y,ax,ay);
        end

        function plotSingularities(obj)
            isSingP1Disc = obj.computeIsSingularP1Discontinous();
            connec = obj.meshDisc.connec;
            xDisc  = obj.meshDisc.coord(:,1);
            yDisc  = obj.meshDisc.coord(:,2);
            figure()
            trisurf(connec,xDisc,yDisc,isSingP1Disc)
            view(0,90)
            colorbar
        end

        function isS = computeIsSingularP1Discontinous(obj)
            isSingP0     = obj.isElemSingular;
            isSingP1Disc = obj.mapP0ToP1Discontinous(isSingP0);
            isS = isSingP1Disc;
        end

        function fP1 = mapP0ToP1Discontinous(obj,f)
            nnodeElem = obj.meshDisc.nnodeElem;
            fRepeted = zeros(size(f,1),nnodeElem);
            for iNode = 1:nnodeElem
                fRepeted(:,iNode) = f;
            end
            fRepeted = transpose(fRepeted);
            fP1 = fRepeted(:);
        end

        function fP1 = mapP1ToP1Discontinous(obj,f)
            nnodeElem = obj.mesh.nnodeElem;
            nElem     = obj.mesh.nelem;
            nDim      = size(f,2);
            fP1 = zeros(nElem,nDim,nnodeElem);
            for iNode = 1:nnodeElem
                nodeI = obj.mesh.connec(:,iNode);
                fP1(:,:,iNode) = f(nodeI,:);
            end
        end        

    end

    methods (Access = private, Static)

        function ab = scalarProduct(a,b)
            ab = a(:,1).*b(:,1) + a(:,2).*b(:,2);
        end

    end

end