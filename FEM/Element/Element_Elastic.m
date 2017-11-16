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
                [obj.B, Bmat] = obj.B.computeB(nunkn,nelem,geometry.nnode,geometry.cartDeriv(:,:,:,igauss));
                
                % Compute Ke
               for iv = 1:geometry.nnode*nunkn
                    for jv = 1:geometry.nnode*nunkn
                        for istre=1:nstre
                            for jstre=1:nstre
                                v = squeeze(Bmat(istre,iv,:).*Cmat(istre,jstre,:).*Bmat(jstre,jv,:));
                                Ke(iv,jv,:) = squeeze(Ke(iv,jv,:)) + v(:).*geometry.area(:,igauss);
                            end
                        end
                    end
                end
            end
            obj.LHS = Ke;
        end
        
        function obj = computeRHS(obj,nunkn,nelem,nnode,bc,idx)
            Fext = zeros(nnode*nunkn,1,nelem);
            for i = 1:length(bc.iN)
                for j = 1:nelem
                    ind = find(idx(:,j) == bc.iN(i));
                    if ~isempty(ind)
                        Fext(ind,:,j) = bc.neunodes(i,3);
                    end
                    % clear ind
                    ind = [];
                end
            end
            obj.RHS = Fext;
        end
        
    end
    
end
