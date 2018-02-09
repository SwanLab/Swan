classdef Element_Thermal < Element
    %Element_Elastic Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! CONSIDER TO IMPLEMENT A CONSTRUCTOR THAT DEFINES B & C DIMENS AT
    % THE PRE-PROCESS !!
    
    properties
    end
    
    methods (Access = ?Physical_Problem)
        function obj = computeLHS(obj)
            Ke = zeros(obj.nunkn*obj.nnode,obj.nunkn*obj.nnode,obj.nelem);
            % Elastic matrix
            for igauss = 1 :obj.geometry.ngaus
                % Strain-displacement matrix
                [obj.B, Bmat] = obj.B.computeB(obj.nunkn,obj.nelem,obj.nnode,obj.geometry.cartd(:,:,:,igauss));
                
                % Compute Ke
                if obj.nelem < 1000 %Just to reduce test.m compute time TO BE REMOVED
                    for i = 1:obj.nelem
                        Ke(:,:,i) = Ke(:,:,i)+Bmat(:,:,i)'*...
                            Bmat(:,:,i)*obj.geometry.dvolu(i,igauss);
                    end
                else
                    for iv=1:obj.nnode*obj.nunkn
                        for jv=1:obj.nnode*obj.nunkn
                            for istre=1:obj.nstre
                                % for jstre=1:nstre
                                v = squeeze(Bmat(istre,iv,:).*Bmat(istre,jv,:));
                                Ke(iv,jv,:) = squeeze(Ke(iv,jv,:)) + v(:).*obj.geometry.dvolu(:,igauss);
                                %end
                            end
                        end
                    end
                end
            end
            obj.LHS = Ke;
        end
    end
    
    methods (Access = protected)
        function Fext = computePuntualRHS(obj)
            Fext = zeros(obj.nnode*obj.nunkn,1,obj.nelem);
            for i = 1:length(obj.bc.iN)
                for j = 1:obj.nelem
                    ind = find(obj.dof.idx(:,j) == obj.bc.iN(i));
                    if ~isempty(ind)
                        Fext(ind,:,j) = obj.bc.neunodes(i,3);
                    end
                    % clear ind
                    ind = [];
                end
            end
        end
        function Fext = computeSuperficialRHS(obj) %To be donne
            Fext = zeros(obj.nnode*obj.nunkn,1,obj.nelem);
        end
        function Fext = computeVolumetricRHS(obj)%To be done
            Fext = zeros(obj.nnode*obj.nunkn,1,obj.nelem);
            
        end
    end
    
end


