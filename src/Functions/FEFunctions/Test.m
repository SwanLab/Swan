classdef Test < BaseFunction

    properties (Access = private)
        uFun
        iDof
    end
    
    methods (Access = public)
        
        function obj = Test(uFun,i)
            obj.uFun  = uFun;
            obj.mesh  = uFun.mesh;
            obj.ndimf = uFun.ndimf;
            obj.ndimfTotal = uFun.ndimfTotal;
            obj.iDof  = i; 
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

        function test = updateMesh(obj,m)
            uF   = LagrangianFunction.create(m,obj.ndimf,obj.uFun.order);
            test = Test(uF,obj.iDof);
        end

    end

    methods (Access = protected)

        function Ni = evaluateNew(obj,xV)
            u     = obj.uFun;
            ndimf = obj.ndimfTotal;
            node = ceil(obj.iDof/ndimf);
            dim  = obj.iDof - (node-1)*ndimf;
            if iscell(xV) % Sample of xV, not necessary all elements
                xV = cell2mat(xV);
                N  = u.computeShapeFunctions(xV);
                nGauss = size(xV,2);
                nEval = size(xV,3);
                Ni = zeros(ndimf,nGauss,nEval);
                Ni(dim,:,:) = N(node,:,:);
            else % All elements              If statement unified if xV always of ndims=3
                N = u.computeShapeFunctions(xV);
                nGauss = size(xV,2);
                nEval = u.mesh.nelem;
                Ni = zeros(ndimf,nGauss,nEval);
                Ni(dim,:,:) = repmat(N(node,:),[1 1 nEval]);
            end
            % transposedDims = length(obj.ndimf):-1:1;
            % extraDims      = length(obj.ndimf)+(1:2);
            % Ni = permute(reshape(Ni,[u.ndimf,nGauss,nEval]),[transposedDims extraDims]);
            Ni = reshape(Ni,[u.ndimf,nGauss,nEval]);
        end

    end


end