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
        
        function obj = computeRHS(obj,nunkn,nelem,nnode,material,bc,idx,geometry,vstrain)
            computeRHS@Element(obj,nunkn,nelem,nnode,bc,idx);
            RHSStrain = obj.computeStrainRHS(nunkn,nelem,material,bc,idx,geometry,vstrain);
            obj.RHS = obj.RHS + RHSStrain;
        end
    end
    
    methods (Access = private)
        function F = computeStrainRHS(obj,nunkn,nelem,material,bc,idx,geometry,vstrain)
            nstre=3;
            Cmat = material.C;
            eforce = zeros(nunkn*geometry.nnode,nunkn*geometry.nnode,nelem);
            stre0=zeros(nstre,1,nelem);
            for igaus=1:geometry.ngaus
                [obj.B, Bmat] = obj.B.computeB(nunkn,nelem,geometry.nnode,geometry.cartDeriv(:,:,:,igaus));
                for iv=1:geometry.nnode*nunkn
                    for jv=1:geometry.nnode*nunkn                        
                        for istre=1:nstre
                            for jstre=1:nstre
%                                  stre0(istre,:) = squeeze(stre0(istre,iv,:)) + squeeze(Cmat(istre,jstre,:)*vstrain(jstre))';
 %                                 eforce(iv,:)=eforce(iv,:)+Bmat(istre,iv,:).*stre0(istre,:).*geometry.dvolu(:,igaus);
                            end
                        end
                    end
                end
            end
            F = eforce;
        end
    end
end
