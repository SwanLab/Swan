classdef Voigt < handle
    
    properties (Access = private)
        A
    end
    
    methods (Access = public)
        
        function obj = Voigt(A)
            obj.init(A)
        end
        
        function voigtA = evaluate(obj,xV)
            matA = obj.A.evaluate(xV);
            ndim = size(matA,1);
            switch ndim
                case 2
                    voigtA = obj.applyVoigt2D(matA);
                case 3
                    voigtA = obj.applyVoigt3D(matA);
            end
        end
    end
    
    methods (Access = private)
        
        function init(obj,A)
            obj.A = A;
        end
        
        function voigtA = applyVoigt2D(obj, matA)
            nPoints = size(matA,3);
            nElem = size(matA,4);
            voigtA = zeros(3,nPoints,nElem);
            
            voigtA(1,:,:) = matA(1,1,:,:); % xx
            voigtA(2,:,:) = matA(2,2,:,:); % yy
            voigtA(3,:,:) = matA(1,2,:,:) + matA(2,1,:,:); % xy

        end
        
        
        
        function voigtA = applyVoigt3D(obj, matA)
            nPoints = size(matA,3);
            nElem = size(matA,4);
            voigtA = zeros(6,nPoints,nElem);
            
            voigtA(1,:,:) = matA(1,1,:,:); % xx
            voigtA(2,:,:) = matA(2,2,:,:); % yy
            voigtA(3,:,:) = matA(3,3,:,:); % zz
            voigtA(4,:,:) = matA(1,2,:,:) + matA(2,1,:,:); % xy
            voigtA(5,:,:) = matA(1,3,:,:) + matA(3,1,:,:); % xz
            voigtA(6,:,:) = matA(2,3,:,:) + matA(3,2,:,:); % yz
        end
        
    end
end
    
    
