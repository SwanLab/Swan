classdef DDP < handle
    % Double Dot Product
    % C = Ai:Bi
    % Ci = Aij:Bj
    properties (Access = private)
        A
        B
    end

    methods (Access = public)

        function obj = DDP(A,B)
            obj.init(A,B)
        end

        function C = evaluate(obj,xV)
            tensorA = obj.A.evaluate(xV);
            tensorB = obj.B.evaluate(xV);
            ndimA = length(size(tensorA));
            ndimB = length(size(tensorB));
            if ndimA == ndimB
                C = obj.computeVectorDDP(tensorA,tensorB);
            elseif ndimA > ndimB
                C = obj.computeMatrixVectorDDP(tensorA,tensorB);
            else
                C = obj.computeMatrixVectorDDP(tensorB,tensorA);
            end

        end

    end

    methods (Access = private)

        function init(obj,A,B)
            obj.A = A;
            obj.B = B;
        end

        function C = computeVectorDDP(obj,A,B)
            C = zeros([1,size(A,[2 3])]);
            nstre = size(A,1);
            for iStre = 1:nstre
                Ai = squeeze(A(iStre,:,:));
                Bi= squeeze(B(iStre,:,:));
                C(1,:,:) = C(1,:,:) + Ai.*Bi;
            end
        end

        function C = computeMatrixVectorDDP(obj,A,B)
            C = zeros(size(A,[1 3 4]));
            nstre = size(A,1);
            for iStre = 1:nstre
                for jStre=1:nstre
                    Aij = squeezeParticular(A(iStre,jStre,:,:),[1 2]);
                    Bj= squeezeParticular(B(jStre,:,:),1);
                    C(iStre,:,:) = squeezeParticular(C(iStre,:,:),1) + Aij.*Bj;
                end
            end
        end

    end

end