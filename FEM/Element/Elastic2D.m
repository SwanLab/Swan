classdef Elastic2D < ElasticDim
    
    properties (Access = protected)
        nstre = 3;
    end
    
    methods (Access = public)
       
        function [B] = computeB(obj,igaus)
            B = zeros(obj.nstre,obj.nnode*obj.dof.nunkn,obj.nelem);
            for i = 1:obj.nnode
                j = obj.dof.nunkn*(i-1)+1;
                B(1,j,:)  = obj.geometry.cartd(1,i,:,igaus);
                B(2,j+1,:)= obj.geometry.cartd(2,i,:,igaus);
                B(3,j,:)  = obj.geometry.cartd(2,i,:,igaus);
                B(3,j+1,:)= obj.geometry.cartd(1,i,:,igaus);
            end
        end        
        
        
    end
    
    methods (Access = protected)
 

                
        function strain = computeStrain(obj,u,idx)
            strain = obj.computeStrain@ElasticDim(u,idx);
            %strain = obj.computeEz(strain,obj.nstre,obj.nelem,obj.material);
        end
               
    end        
    
    methods (Access = private, Static)
        
        function strain = computeEz(strain,nstre,nelem,material)
            mu = material.mu;
            kappa = material.kappa;
            epoiss = (kappa(1,1) - mu(1,1))./(kappa(1,1) + mu(1,1));
            epoiss = full(ones(1,nelem)*epoiss);
            strain(nstre+1,:,:) = (-epoiss./(1-epoiss)).*(strain(1,:,:)+strain(2,:,:));
        end
        
    end
end

