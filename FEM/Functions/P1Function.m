classdef P1Function < FeFunction

    properties (Access = public)
    end

    properties (Access = private)
        interpolation
    end

    properties (Access = private)
        type
        connec
    end

    methods (Access = public)

        function obj = P1Function(cParams)
            obj.init(cParams);
            obj.createInterpolation();
        end

        function fxV = evaluate(obj, xV)
            obj.interpolation.computeShapeDeriv(xV);
            shapes = obj.interpolation.shape;
            nNode  = size(shapes,1);
            nGaus  = size(shapes,2);
            nF     = size(obj.fValues,2);
            nElem  = size(obj.connec,1);
            fxV = zeros(nF,nGaus,nElem);
            for iGaus = 1:nGaus
                for iNode = 1:nNode
                    node = obj.connec(:,iNode);
                    Ni = shapes(iNode,iGaus,:);
                    fi = obj.fValues(node,:);
                    f(:,1,:) = Ni*fi';
                    fxV(:,iGaus,:) = fxV(:,iGaus,:) + f;
                end
            end

        end

        function plot(obj, m) % 2D domains only
            dim = 1;
            x = m.coord(:,1);
            y = m.coord(:,2);
            z = obj.fValues(:,dim);
            figure()
            a = trisurf(m.connec,x,y,z);
            view(0,90)
%             colorbar
            shading interp
            a.EdgeColor = [0 0 0];
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.connec  = cParams.connec;
            obj.type    = cParams.type;
            obj.fValues = cParams.fValues;
            obj.ndimf   = size(cParams.fValues,2);
        end

        function createInterpolation(obj)
            m.type = obj.type;
            obj.interpolation = Interpolation.create(m,'LINEAR');
        end

        function fRep = repeatFunctionAtNodes(obj)
           f = obj.fValues;
           nNode  = size(obj.connec,2);
           nDime  = size(f,2);
           nElem  = size(obj.connec,1);
           fNodeElem = zeros(nDime,nNode,nElem);
           fNods  = transpose(f);
           for inode = 1:nNode
               nodes = obj.connec(:,inode);
               fNode = fNods(:,nodes);
               fNodeElem(:,inode,:) = fNode;
           end
           fRep = fNodeElem;
        end

    end

end