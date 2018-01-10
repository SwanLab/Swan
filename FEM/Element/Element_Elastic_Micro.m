classdef Element_Elastic_Micro < Element_Elastic
    %Element_Elastic Summary of this class goes here
    %   Detailed explanation goes here
    % !! CONSIDER TO IMPLEMENT A CONSTRUCTOR THAT DEFINES B & C DIMENS AT
    % THE PRE-PROCESS !!
    
    properties
    end
    
    methods (Access = {?Physical_Problem, ?Element})
        function obj = Element_Elastic_Micro
            obj.B = B2;
        end
        
        function obj = computeRHS(obj,nunkn,nstre,nelem,nnode,material,bc,idx,geometry,vstrain)
            computeRHS@Element(obj,nunkn,nelem,nnode,bc,idx);
            RHSStrain = obj.computeStrainRHS(nunkn,nstre,nelem,material,bc,idx,geometry,vstrain);
            obj.RHS = obj.RHS + RHSStrain;             
        end
    end
    
    methods (Access = private)
        function F = computeStrainRHS(obj,nunkn,nstre,nelem,material,bc,idx,geometry,vstrain)
            Cmat = material.C;            
            eforce = zeros(nunkn*geometry.nnode,1,nelem);
            sigma=zeros(nstre,1,nelem);
            for igaus=1:geometry.ngaus
                [obj.B, Bmat] = obj.B.computeB(nunkn,nelem,geometry.nnode,geometry.cartDeriv(:,:,:,igaus));
                for istre=1:nstre
                    for jstre=1:nstre
                        sigma(istre,:) = sigma(istre,:) + squeeze(Cmat(istre,jstre,:)*vstrain(jstre))';
                    end
                end                
                for iv=1:geometry.nnode*nunkn
                    for istre=1:nstre
                        eforce(iv,:)=eforce(iv,:)+(squeeze(Bmat(istre,iv,:)).*sigma(istre,:)'.*geometry.dvolu(:,igaus))';
                    end
                end
            end
            F = -eforce;
        end
    end
end
