classdef FeFunction < handle
    
   properties (Access = public)
       fElem
   end
    
   properties (Access = private)
       interpolation
   end
   
   properties (Access = private)
      mesh 
      fNodes
   end
    
   methods (Access = public)
      
       function obj = FeFunction(cParams)
           obj.init(cParams);
           obj.createInterpolation();
           obj.computeFnodesByelem();
       end
       
       function fC = computeValueInCenterElement(obj)          
            m = obj.mesh;
            q = Quadrature.set(m.geometryType);
            q.computeQuadrature('CONSTANT');
            xV = q.posgp;
            fCenter = obj.interpolateFunction(xV);
            fC = squeeze(fCenter);            
       end
        
        function computeFnodesByelem(obj)
            f = obj.fNodes;
            nNode  = obj.mesh.nnode;
            nDime  = size(f,2);
            nElem  = obj.mesh.nelem;
            coordE = zeros(nDime,nNode,nElem);
            coords = transpose(f);
            for inode = 1:nNode
                nodes = obj.mesh.connec(:,inode);
                coordNodes = coords(:,nodes);
                coordE(:,inode,:) = coordNodes;
            end
            obj.fElem = coordE;           
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
       
   end
   
   methods (Access = private)
      
       function init(obj,cParams)
           obj.mesh   = cParams.mesh;
           obj.fNodes = cParams.fNodes;           
       end
       
       function createInterpolation(obj)
           obj.interpolation = Interpolation.create(obj.mesh,'LINEAR');
       end
       

   end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end