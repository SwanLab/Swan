classdef Element_Thermal < Element
    %Element_Elastic Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! CONSIDER TO IMPLEMENT A CONSTRUCTOR THAT DEFINES B & C DIMENS AT
    % THE PRE-PROCESS !!
    
    properties
    end
    
    methods (Access = ?Physical_Problem)
        function [r,dr] = computeResidual(obj,uL)
            % *************************************************************
            % Compute
            % - residual: r = Ku - F
            % - residual derivative: dr = K
            % *************************************************************
            [K] = obj.computeStiffnessMatrix();
                      
            % Assemble
            [K] = obj.AssembleMatrix(K);
            
            %Assemble u and Fext
            u = zeros(obj.dof.ndof,1);
            u(obj.dof.vL) = uL;
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
        
        function [K] = computeStiffnessMatrix(obj)
            
            % Stiffness matrix
            Ke = zeros(obj.nunkn*obj.nnode,obj.nunkn*obj.nnode,obj.nelem);
            
            obj.B.value = cell(obj.geometry.ngaus);
            for igaus = 1 :obj.geometry.ngaus
                % Strain-displacement matrix
                Bmat = obj.B.computeB(obj.nunkn,obj.nelem,obj.nnode,obj.geometry.cartd(:,:,:,igaus));
                
                % Compute Ke
                if obj.nelem < 1000 %Just to reduce test.m compute time TO BE REMOVED
                    for i = 1:obj.nelem
                        Ke(:,:,i) = Ke(:,:,i)+Bmat(:,:,i)'*Bmat(:,:,i)*obj.geometry.dvolu(i,igaus);
                    end
                else
                    for iv = 1:obj.nnode*obj.nunkn
                        for jv = 1:obj.nnode*obj.nunkn
                            for istre = 1:obj.nstre
                                % for jstre=1:nstre
                                v = squeeze(Bmat(istre,iv,:).*Bmat(istre,jv,:));
                                Ke(iv,jv,:) = squeeze(Ke(iv,jv,:)) + v(:).*obj.geometry.dvolu(:,igaus);
                                %end
                            end
                        end
                    end
                end
                obj.B.value{igaus} = Bmat;
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


