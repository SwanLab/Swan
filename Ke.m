classdef Ke<Matrix_Elemental
    %UNTITLED9 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = Ke(nstre,nunkn,nelem,geometry)
            obj.value=zeros(nunkn*geometry.nnode,nunkn*geometry.nnode,nelem);
            % Elastic matrix
            Cmat = C(nstre,nelem,geometry.ndime);
            
            for igauss=1:geometry.ngaus
            % Strain-displacement matrix
            Bmat = B(nstre,nunkn,nelem,geometry.nnode,geometry.cartDeriv(:,:,:,igauss),geometry.ndime);
            
            % Compute Ke
            for i = 1:nelem
               obj.value(:,:,i) = obj.value(:,:,i)+Bmat.value(:,:,i)'*Cmat.value(:,:,i)*...
                   Bmat.value(:,:,i)*geometry.area(i,igauss);
            end
            end

        end
    end
    
end

