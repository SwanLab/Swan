classdef Element_Elastic < Element
    %Element_Elastic Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! CONSIDER TO IMPLEMENT A CONSTRUCTOR THAT DEFINES B & C DIMENS AT
    % THE PRE-PROCESS !!
    
    properties
    end
    
    methods (Access = ?Physical_Problem)
        function obj = computeLHS(obj,nunkn,nstre,nelem,geometry,material)
            
            Ke = zeros(nunkn*geometry.nnode,nunkn*geometry.nnode,nelem);
            % Elastic matrix
             Cmat = material.C;           
            
            for igauss = 1 :geometry.ngaus
                % Strain-displacement matrix
                [obj.B, Bmat] = obj.B.computeB(nunkn,nelem,geometry.nnode,geometry.cartd(:,:,:,igauss));
                
                for iv=1:geometry.nnode*nunkn
                    for jv=1:geometry.nnode*nunkn
                        for istre=1:nstre
                            for jstre=1:nstre
                                v = squeeze(Bmat(istre,iv,:).*Cmat(istre,jstre,:).*Bmat(jstre,jv,:));
                                Ke(iv,jv,:) = squeeze(Ke(iv,jv,:)) + v(:).*geometry.dvolu(:,igauss);
                            end
                        end
                        
                    end
                end
            end
            obj.LHS = Ke;
        end
    end
    methods (Access = protected)
        function Fext = computeSuperficialRHS(obj,nunkn,nelem,nnode,bc,idx) %To be donne
            Fext = zeros(nnode*nunkn,1,nelem);
        end
        function Fext = computeVolumetricRHS(obj,nunkn,nelem,nnode,bc,idx)%To be done
            Fext = zeros(nnode*nunkn,1,nelem);
            
        end
    end
    
end
