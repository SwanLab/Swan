classdef Test < handle
    
    properties (Access = public)
                mesh
            ndimf
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        uFun
        iDof
    end
    
    methods (Access = public)
        
        function obj = Test(uFun,i)
            obj.uFun = uFun;
            obj.mesh = uFun.mesh;
            obj.ndimf = uFun.ndimf;
            obj.iDof = i; 
        end
        
        function gradN = computeGrad(obj,xV)
            u = obj.uFun;
            dNdx = u.evaluateCartesianDerivatives(xV);
            ndim = u.mesh.ndim;
            node = ceil(obj.iDof/ndim);
            dim  = obj.iDof - (node-1)*ndim;
            
            nGauss = size(xV,2);
            nElem  = u.mesh.nelem;
            gradN = zeros(ndim,ndim,nGauss,nElem);
            gradN(dim,:,:,:) = dNdx(:,node,:,:);
   
        end
        
    end
  
    
end