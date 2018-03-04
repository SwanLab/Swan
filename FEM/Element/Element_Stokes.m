classdef Element_Stokes < Element
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        LHS
        RHS
    end
    
    methods (Access = ?Physical_Problem)
        
        function obj = compute_LHS(obj,dim,nelem,geometry_variable,material,nfields)
            for ifield = 1:nfields
                for jfields = 1:nfields
                    obj.LHS{ifield,jfield} = obj.computeMatrix(dim,nelem,geometry_variable(ifield),geometry_variable(jfield),material,ifield,jfield);
                end
            end
        end
        
        function obj = compute_RHS(obj,dim,nelem,geometry_variable,dof)
            
        end

        function K = compute_K(obj,dim,nelem,geometry_test,geometry_variable,material)
             nunkn = dim.nunkn(1);
            K = zeros(nunkn*geometry_test.nnode,nunkn*geometry_variable.nnode,nelem);
           
            Cmat = material.mu;
            
            for igauss = 1 :geometry_variable.ngaus
                
                Bmat = obj.computeB(nunkn,nelem,geometry_variable.nnode,geometry_variable.cartDeriv(:,:,:,igauss));
                
%                 B_p=reshape(Bmat,[geometry.nnode*nunkn,1,nelem]);
                
                for iv=1:geometry_test.nnode*nunkn  
                    for jv=1:geometry_variable.nnode*nunkn
                        for istre=1:nunkn*geometry_test.ndime
                            for jstre=1:nunkn*geometry_variable.ndime
                                v = squeeze(Bmat(istre,iv,:).*Cmat(istre,jstre,:).*Bmat(jstre,jv,:));
                                K(iv,jv,:)= squeeze(K(iv,jv,:)) + v(:).*geometry_variable.dvolu(:,igauss);
%                                 obj.LHS(iv,jv,:) = squeeze(obj.LHS(iv,jv,:)) + v(:).*geometry.dvolu(:,igauss);
                            end
                        end
                    end
                end
                
            end
            
          
            
            end
        
        function D = compute_D(obj,dim,nelem,geometry_test,geometry_variable)
            nunkn_u=dim.nunkn(1);
 
            
            D = zeros(nunkn_u*geometry_test.nnode,geometry_variable.nnode,nelem);

            for igauss=1:geometry_test.ngaus
                for inode_var = 1:geometry_variable.nnode
                    for inode_test = 1:geometry_test.nnode
                        for idime = 1:geometry_test.ndime
                            dof_test = inode_test*nunkn_u - nunkn_u + idime;
                            v= squeeze (geometry_test.cartDeriv(idime,inode_test,:,igauss));
                            D(dof_test,inode_var,:)= squeeze(D(dof_test,inode_var,:)) - v(:).*geometry_variable.shape(inode_var,igauss)...
                                .*geometry_test.dvolu(:,igauss);
                            
                        end
                    end
                end
            end
        end
        
       function B = computeB(obj,nunkn,nelem,nnode,cartd)
            B = zeros(2,nnode*nunkn,nelem);
            for i = 1:nnode
                j = nunkn*(i-1)+1;
                B(1,j,:)  = cartd(1,i,:);
                B(2,j+1,:)= cartd(1,i,:);
                B(3,j,:)  = cartd(2,i,:);
                B(4,j+1,:)= cartd(2,i,:);
            end
        end
        
        
    end
    methods (Access=protected)
        function mat = computeMatrix(obj,dim,nelem,geometry_test,geometry_variable,material,ifield,jfield)

           if ifield == 1 && jfield==1
               mat = obj.compute_K(dim,nelem,geometry_test,geometry_variable,material);
           elseif ifield == 1 && jfield==2
               mat = obj.compute_D(dim,nelem,geometry_test,geometry_variable);
           elseif ifield == 2 && jfield==1
               mat = permute(obj.LHS{1,2},[2,1,3]);
           else
%                D = obj.LHS{2,1};
%                D_traspose= obj.LHS{1,2};
%                row = length(D(:,1,1));
%                col = length (D_traspose(1,:,1));
               D = obj.LHS{1,1};
               D_traspose= obj.LHS{1,1};
               row = length(D(:,1,1));
               col = length (D_traspose(1,:,1));

               mat = zeros(row,col,nelem);
           end         

%             obj.LHS= [[K D]; [permute(D,[2,1,3]) zeros(1,1,nelem)]];
        end
        
        function Fext = computePuntualRHS(obj,nunkn,nelem,nnode,bc,idx)
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
        end
        function Fext = computeSuperficialRHS(obj,nunkn,nelem,nnode,bc,idx) %To be donne
%             Fext = zeros(nnode*nunkn,1,nelem);
              Fext=0;
        end
        function Fext = computeVolumetricRHS(obj,dim,nelem,geometry_variable,bc,dof)%To be done
            idx = dof(1).idx;
             geometry_variable = geometry_variable(1);
            nnode = geometry_variable.nnode;
            nunkn= dim.nunkn(1);
%             f = zeros(nnode*nunkn,1,nelem);
            Fext = zeros(nnode*nunkn,1,nelem);
            obj.RHS = zeros(nnode*nunkn,1,nelem);
%             for i = 1:length(bc.iN)
%                 for j = 1:nnode*nunkn
%                     ind = find(idx(j,:) == bc.iN(i));
%                     if ~isempty(ind)
%                         f(j,:,ind) = bc.force(i,3);
%                     end
%                     %                     clear ind
%                     ind = [];
%                 end
%             end

             
%              for igaus=1:geometry_variable.ngaus
%                 for inode=1:nnode
%                     for jnode=1:nnode
%                         for iunkn=1:nunkn
%                             elemental_dof = jnode*nunkn-nunkn+iunkn; %% dof per guardar el valor de la integral
% 
%                                 v= squeeze(geometry_variable.shape(inode,igaus).*geometry_variable.shape(jnode,igaus).*f(elemental_dof,1,:));
%                                 Fext(elemental_dof,1,:)= squeeze(Fext(elemental_dof,1,:)) + v(:).*geometry_variable.dvolu(:,igaus);
%                             
%                         end
%                     end
%                 end
%             end

            if  ~isempty(bc.force)
                f=bc.force;

                for igaus=1:geometry_variable.ngaus
                    for inode=1:nnode
                        for iunkn=1:nunkn
                            elemental_dof = inode*nunkn-nunkn+iunkn; %% dof per guardar el valor de la integral

                            v= squeeze(geometry_variable.shape(inode,igaus).*f(iunkn,igaus,:));
                            Fext(elemental_dof,1,:)= squeeze(Fext(elemental_dof,1,:)) + v(:).*geometry_variable.dvolu(:,igaus);

                        end
                    end
                end
            end

        end
    end
    
end