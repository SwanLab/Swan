classdef Element_Elastic < Element
    %Element_Elastic Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! CONSIDER TO IMPLEMENT A CONSTRUCTOR THAT DEFINES B & C DIMENS AT
    % THE PRE-PROCESS !!
    
    properties
    end
    
    methods (Access = {?Physical_Problem, ?Element})
        function obj = Element_Elastic()
            obj.nincr = 1;
            obj.cload = 0;
        end
        
        function [r,dr] = computeResidual(obj,uL)
            % *************************************************************
            % Compute
            % - residual: r = Ku - F
            % - residual derivative: dr = K
            % *************************************************************
            % Compute stiffness matrix
            K = obj.computeStiffnessMatrix();
                     
            % Assemble
            K = obj.AssembleMatrix(K);
            
            % Assemble u and Fext
            u = zeros(obj.dof.ndof,1);
            u(obj.dof.vL) = uL;
            if ~isempty(obj.dof.vR)
                u(obj.dof.vR) = obj.bc.fixnodes(:,3);
                fext = obj.cload(obj.dof.vL)-K(obj.dof.vL,obj.dof.vR)*u(obj.dof.vR); % fext + reac
            else
                fext = obj.cload(obj.dof.vL);
            end
            fint = K(obj.dof.vL,obj.dof.vL)*u(obj.dof.vL);
            r = fint - fext;
            dr = K(obj.dof.vL, obj.dof.vL);
        end
                
        function [K] = computeStiffnessMatrix(obj)
            
            % Stiffness matrix
            Ke = zeros(obj.nunkn*obj.nnode,obj.nunkn*obj.nnode,obj.nelem);
            
            % Elastic matrix
            Cmat = obj.material.C;
            
            obj.B.value = cell(obj.geometry.ngaus);
            for igaus = 1 :obj.geometry.ngaus
                % Strain-displacement matrix
                Bmat = obj.B.computeB(obj.nunkn,obj.nelem,obj.nnode,obj.geometry.cartd(:,:,:,igaus));
                
                for iv = 1:obj.nnode*obj.nunkn
                    for jv = 1:obj.nnode*obj.nunkn
                        for istre = 1:obj.nstre
                            for jstre = 1:obj.nstre
                                v = squeeze(Bmat(istre,iv,:).*Cmat(istre,jstre,:).*Bmat(jstre,jv,:));
                                Ke(iv,jv,:) = squeeze(Ke(iv,jv,:)) + v(:).*obj.geometry.dvolu(:,igaus);
                            end
                        end
                        
                    end
                end
                obj.B.value{igaus} = Bmat;
            end
            K = Ke;
        end        
        
    end
    
    
    methods(Access = protected) % Only the child sees the function
        function FextSuperficial = computeSuperficialFext(obj,bc)
            FextSuperficial = zeros(obj.nnode*obj.nunkn,1,obj.nelem);
        end
        
        function FextVolumetric = computeVolumetricFext(obj,bc)
            FextVolumetric = zeros(obj.nnode*obj.nunkn,1,obj.nelem);
        end
        
        function variables = computeDispStressStrain(obj,uL)
            variables.d_u = zeros(obj.dof.ndof,1);
            variables.d_u(obj.dof.vL) = uL;
            variables.d_u(obj.dof.vR) = obj.bc.fixnodes(:,3);
            variables.strain = obj.computeStrain(variables.d_u,obj.dim,obj.nnode,obj.nelem,obj.geometry.ngaus,obj.dof.idx);
            variables.stress = obj.computeStress(variables.strain,obj.material.C,obj.geometry.ngaus,obj.nstre);
        end
        
          function strain = computeStrain(obj,d_u,dim,nnode,nelem,ngaus,idx)
            strain = zeros(dim.nstre,nelem,ngaus);
            for igaus = 1:ngaus
                Bmat = obj.B.value{igaus};
                %Bmat = Bmat{1,1};
                for istre=1:dim.nstre
                    for inode=1:nnode
                        for idime=1:dim.nunkn
                            ievab = dim.nunkn*(inode-1)+idime;
                            strain(istre,:,igaus)=strain(istre,:,igaus)+(squeeze(Bmat(istre,ievab,:)).*d_u(idx(ievab,:)))';
                        end
                    end
                end
            end
        end

    end
    
    methods(Static)
        function variables = permuteStressStrain(variables)
            variables.strain = permute(variables.strain, [3 1 2]);
            variables.stress = permute(variables.stress, [3 1 2]);
        end
    end
    
    methods(Static, Access = protected)
        % Only used in Element_Elastic_2D
        function strain = computeEz(strain,nstre,nelem,material)
            mu = material.mu;
            kappa = material.kappa;
            epoiss = (kappa(1,1) - mu(1,1))./(kappa(1,1) + mu(1,1));
            epoiss = ones(1,nelem)*epoiss;
            strain(nstre+1,:,:) = (-epoiss./(1-epoiss)).*(strain(1,:,:)+strain(2,:,:));
        end
        
        % Compute strains (e = B�u)
      
        
        % Compute stresses
        function stres = computeStress(strain,C,ngaus,nstre)
            stres = zeros(size(strain));
            for igaus = 1:ngaus
                for istre=1:nstre
                    for jstre=1:nstre
                        stres(istre,:,igaus) = stres(istre,:,igaus) + squeeze(C(istre,jstre,:))'.*strain(jstre,:,igaus);
                    end
                end
            end
        end
    end
end
