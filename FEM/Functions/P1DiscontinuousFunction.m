classdef P1DiscontinuousFunction < FeFunction
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        connec
        fNodes
    end
    
    properties (Access = private)
    end
    
    methods (Access = public)
        
        function obj = P1DiscontinuousFunction(cParams)
            obj.init(cParams)
            obj.computeFvalues();
        end

        function fV = evaluate(obj, xV)
            % Goal: develop this function
            func = obj.fValues;
            obj.interpolation.computeShapeDeriv(xV);
            shapes = obj.interpolation.shape;
            nNode  = size(shapes,1);
            nGaus  = size(shapes,2);
            nF     = size(func,1);
            nElem  = size(func,3);
            fxV = zeros(nF,nGaus,nElem);
            for kNode = 1:nNode
                shapeKJ = shapes(kNode,:,:);
                fKJ     = func(:,kNode,:);
                f = bsxfun(@times,shapeKJ,fKJ);
                fxV = fxV + f;
            end
        end

        function plot(obj, m)
            dim = 1;
            ndims   = size(obj.fValues, 1);
            nelem   = size(obj.connec, 1);
            nnodeEl = size(obj.connec, 2);
            mD = m.createDiscontinuousMesh();
            fD = reshape(obj.fValues, [ndims, nelem*nnodeEl])';
            figure()
            trisurf(mD.connec, mD.coord(:,1), mD.coord(:,2), fD(:,dim))
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fNodes = cParams.fNodes;
            obj.connec = cParams.connec;
        end

        function computeFvalues(obj)
           f = obj.fNodes;
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
           obj.fValues = fNodeElem; 
        end
        
    end
    
end