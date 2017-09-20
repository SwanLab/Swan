classdef Ke<Matrix_Elemental
    %UNTITLED9 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = Ke(nstre,nunkn,nelem,geometry)
            
            % Strain-displacement matrix
            Bmat = B(nstre,nunkn,nelem,geometry.nnode,geometry.cartDeriv);
            
            % Elastic matrix
            Cmat = C(nstre,nelem);
            
            % Compute Ke
            for i = 1:nelem
               obj.value(:,:,i) = Bmat.value(:,:,i)'*Cmat.value(:,:,i)*...
                   Bmat.value(:,:,i)*geometry.area(i);
            end

        end
    end
    
end

