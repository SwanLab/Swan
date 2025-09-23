classdef Test < BaseFunction
    
    properties (Access = public)
          %      mesh
          %  ndimf
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
            ndimf = u.ndimf;
            node = ceil(obj.iDof/ndimf);
            dim  = obj.iDof - (node-1)*ndimf;            
            nGauss = size(xV,2);
            nElem  = u.mesh.nelem;
            gradN = zeros(ndimf,u.mesh.ndim,nGauss,nElem);
            gradN(dim,:,:,:) = dNdx(:,node,:,:);
   
        end

    end

    methods (Access = protected)

        function Ni = evaluateNew(obj,xV)
            u     = obj.uFun;
            N  = u.computeShapeFunctions(xV);
            ndimf = u.ndimf;
            node = ceil(obj.iDof/ndimf);
            dim  = obj.iDof - (node-1)*ndimf;
            nGauss = size(xV,2);
            if ismatrix(xV)
                nEval = u.mesh.nelem;
                Ni = zeros(ndimf,nGauss,nEval);
                Ni(dim,:,:) = repmat(N(node,:),[1 1 nEval]);
            else
                nEval = size(xV,3);
                Ni = zeros(ndimf,nGauss,nEval);
                Ni(dim,:,:) = N(node,:,:);
            end
        end

    end


end