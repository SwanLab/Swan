classdef P1Function < FeFunction
    
   properties (Access = public)
       fElem
   end
    
   properties (Access = private)
       interpolation
   end
   
   properties (Access = private)
      type
      connec
      fNodes
   end
    
   methods (Access = public)
      
       function obj = P1Function(cParams)
           obj.init(cParams);
           obj.createInterpolation();
           obj.computeFnodesByelem();
       end
       
       function fC = computeValueInCenterElement(obj)
            q = Quadrature.set(obj.type);
            q.computeQuadrature('CONSTANT');
            xV = q.posgp;
            fCenter = obj.interpolateFunction(xV);
            fC = squeeze(fCenter);
       end
        
       function fxV = interpolateFunction(obj,xV)
            func = obj.fElem;
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

       function fDisc = computeDiscontinousField(obj)
            nodesCont  = obj.connec';
            nodesCont  = nodesCont(:);
            fDisc  = zeros(size(nodesCont));
            f      = obj.fNodes;
            for iDir = 1:size(f,2)
              fCont = f(nodesCont,iDir);
              fDisc(:,iDir) = fCont;
            end
       end
       
       function plot(obj, m) % 2D domains only
           figure()
           trisurf(m.connec, m.coord(:,1), m.coord(:,2), obj.fNodes(:,1))
       end

   end
   
   methods (Access = private)
      
       function init(obj,cParams)
           obj.connec = cParams.connec;
           obj.type   = cParams.type;
           obj.fNodes = cParams.fNodes;
       end
       
       function createInterpolation(obj)
           m.type = obj.type;
           obj.interpolation = Interpolation.create(m,'LINEAR');
       end
       
       function computeFnodesByelem(obj)
           f = obj.fNodes;
           nNode  = size(obj.connec,2);
           nDime  = size(f,2);
           nElem  = size(obj.connec,1);
           coordE = zeros(nDime,nNode,nElem);
           coords = transpose(f);
           for inode = 1:nNode
               nodes = obj.connec(:,inode);
               coordNodes = coords(:,nodes);
               coordE(:,inode,:) = coordNodes;
           end
           obj.fElem = coordE;
       end

   end

end