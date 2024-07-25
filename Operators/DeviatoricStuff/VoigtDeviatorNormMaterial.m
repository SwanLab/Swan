classdef VoigtDeviatorNormMaterial < handle

properties (Access = private)
    nElements
    nStresses
end

    methods (Access = public)
        function obj = VoigtDeviatorNormMaterial(N,nEval)
            obj.nElements = nEval;
            switch N
                case 2
                    obj.nStresses = 3;
                case 3
                    obj.nStresses = 6;
            end
        end

        function A = evaluate(obj,xV)
            nGaus = size(xV,2);
            nElem = obj.nElements;
            nStre = obj.nStresses;
            A     = zeros(nStre,nStre,nGaus,nElem);
            switch nStre
                case 3
                    A(1,1,:,:) = 2;
                    A(2,2,:,:) = 2;
                    A(3,3,:,:) = 1;
                case 6
                    A(1,1,:,:) = 2;
                    A(2,2,:,:) = 2;
                    A(3,3,:,:) = 2;
                    A(4,4,:,:) = 1;
                    A(5,5,:,:) = 1;
                    A(6,6,:,:) = 1;
            end
        end
    end
end