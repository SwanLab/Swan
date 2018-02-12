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
        
        function obj = computeRHS(obj,vstrain)
            computeRHS@Element(obj);
            RHSStrain = obj.computeStrainRHS(vstrain);
            obj.RHS = obj.RHS + RHSStrain;             
        end
    end
    
    methods (Access = private)
        function F = computeStrainRHS(obj,vstrain)
            Cmat = obj.material.C;            
            eforce = zeros(obj.nunkn*obj.nnode,1,obj.nelem);
            sigma=zeros(obj.nstre,1,obj.nelem);
            for igaus=1:obj.geometry.ngaus
                Bmat = obj.B.computeB(obj.nunkn,obj.nelem,obj.nnode,obj.geometry.cartd(:,:,:,igaus));
                for istre=1:obj.nstre
                    for jstre=1:obj.nstre
                        sigma(istre,:) = sigma(istre,:) + squeeze(Cmat(istre,jstre,:)*vstrain(jstre))';
                    end
                end                
                for iv=1:obj.nnode*obj.nunkn
                    for istre=1:obj.nstre
                        eforce(iv,:)=eforce(iv,:)+(squeeze(Bmat(istre,iv,:)).*sigma(istre,:)'.*obj.geometry.dvolu(:,igaus))';
                    end
                end
            end
            F = -eforce;
        end
    end
end
