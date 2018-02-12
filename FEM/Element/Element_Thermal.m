classdef Element_Thermal < Element
    %Element_Elastic Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! CONSIDER TO IMPLEMENT A CONSTRUCTOR THAT DEFINES B & C DIMENS AT
    % THE PRE-PROCESS !!
    
    properties
    end
    
    methods (Access = ?Physical_Problem)
        function [r,dr] = computeResidual(obj,u)
            % *************************************************************
            % Compute
            % - residual: r = Ku - F
            % - residual derivative: dr = K
            % *************************************************************
            [K,B] = obj.computeStiffnessMatrix();
            
            % Stores strain-displacement matrix
            obj.B = B;
            
            % Assemble
            [K] = obj.AssembleMatrix(K);
            
            if ~isempty(obj.dof.vR)
                u(obj.dof.vR) = obj.bc.fixnodes(:,3);
                fext = obj.Fext(obj.dof.vL)-K(obj.dof.vL,obj.dof.vR)*u(obj.dof.vR);
            else
                fext = obj.Fext(obj.dof.vL);
            end
            fint = K(obj.dof.vL,obj.dof.vL)*u(obj.dof.vL);
            r = fint - fext;
            dr = K(obj.dof.vL,obj.dof.vL);
        end
        
        function [K,B] = computeStiffnessMatrix(obj)
            
            % Stiffness matrix
            Ke = zeros(obj.nunkn*obj.nnode,obj.nunkn*obj.nnode,obj.nelem);
            
            for igauss = 1 :obj.geometry.ngaus
                % Strain-displacement matrix
                [B, Bmat] = obj.B.computeB(obj.nunkn,obj.nelem,obj.nnode,obj.geometry.cartd(:,:,:,igauss));
                
                % Compute Ke
                if obj.nelem < 1000 %Just to reduce test.m compute time TO BE REMOVED
                    for i = 1:obj.nelem
                        Ke(:,:,i) = Ke(:,:,i)+Bmat(:,:,i)'*Bmat(:,:,i)*obj.geometry.dvolu(i,igauss);
                    end
                else
                    for iv = 1:obj.nnode*obj.nunkn
                        for jv = 1:obj.nnode*obj.nunkn
                            for istre = 1:obj.nstre
                                % for jstre=1:nstre
                                v = squeeze(Bmat(istre,iv,:).*Bmat(istre,jv,:));
                                Ke(iv,jv,:) = squeeze(Ke(iv,jv,:)) + v(:).*obj.geometry.dvolu(:,igauss);
                                %end
                            end
                        end
                    end
                end
            end
            K = Ke;
        end
    end
    
    methods (Access = protected)
        function FextSuperficial = computeSuperficialFext(obj,bc)
            FextSuperficial = zeros(obj.nnode*obj.nunkn,1,obj.nelem);
        end
        
        function FextVolumetric = computeVolumetricFext(obj,bc)
            FextVolumetric = zeros(obj.nnode*obj.nunkn,1,obj.nelem);
        end
    end
    
end


