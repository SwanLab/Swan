classdef Element_Stokes < Element
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        LHS_elem
        LHS
        RHS
        M_elem
        K_elem
    end
    
    methods (Access = ?Physical_Problem)
        function [r,dr] = computeResidual(obj,x,dr,x_n)
%             K = compute_LHS(obj);
            if (nargin ==3)
                Mred_x_n = zeros(length(obj.dof.free{1}),1);
            else
               M = obj.AssembleMatrix(obj.M_elem,1,1);
               Mred = M(obj.dof.free{1},obj.dof.free{1});
               Mred_x_n = Mred*x_n;
            end
            
            Fext = compute_RHS(obj);
            
%             K = obj.AssembleMatrix(obj.K_elem,1,1);

            R = obj.compute_imposed_displacemet_force(obj.LHS);
            Fext = Fext + R ;
            
%             Kred = obj.full_matrix_2_reduced_matrix(K);

            Fext_red = obj.full_vector_2_reduced_vector(Fext);
            Fext_red(1:length(obj.dof.free{1}),1) = Fext_red(1:length(obj.dof.free{1}),1) + Mred_x_n;
            
            fint_red = dr*x;

            r = fint_red - (Fext_red);
%             dr = Kred;
            
        end
        
        function dr = computedr(obj,dt)
            if nargin < 2
                dt=inf;
            end
             obj.LHS = compute_LHS(obj,dt);
             LHSred = obj.full_matrix_2_reduced_matrix(obj.LHS);
             dr = LHSred;
        end
        
        function LHS = compute_LHS(obj,dt)
            for ifield = 1:obj.nfields
                for jfield = 1:obj.nfields
                    obj.LHS_elem{ifield,jfield} = obj.computeMatrix(obj.dim,obj.nelem,obj.geometry(ifield),obj.geometry(jfield),obj.material,dt,ifield,jfield);
%                         mat = obj.computeMatrix(obj.dim,obj.nelem,obj.geometry(ifield),obj.geometry(jfield),obj.material,ifield,jfield); 
                    LHS{ifield,jfield} = obj.AssembleMatrix(obj.LHS_elem{ifield,jfield},ifield,jfield);
                        
                end
            end
%             LHS= AssembleMatrix(obj,obj.LHS_elem);
            LHS = cell2mat(LHS);
        end
        
        function RHS = compute_RHS(obj)
            Fext = obj.computeVolumetricFext(obj.dim,obj.nelem,obj.geometry,obj.dof);
            g = obj.compute_velocity_divergence;
            RHS_elem{1,1} = Fext;
            RHS_elem{2,1} = g;
            
            RHS = AssembleVector(obj,RHS_elem);
        end

        function M = compute_M(obj,dim,nelem,geometry_test,geometry_variable,dt)
            nunkn = dim.nunkn(1);
            M = zeros(nunkn*geometry_test.nnode,nunkn*geometry_variable.nnode,nelem);
            
             for igauss = 1 :geometry_variable.ngaus       
                for inode= 1:geometry_test.nnode
                    for jnode= 1:geometry_variable.nnode
                        for iunkn= 1:nunkn
                            for junkn= 1:nunkn
                                v = squeeze(geometry_test.shape(inode,igauss,:).*geometry_variable.shape(jnode,igauss,:));
                                M(nunkn*(inode-1)+iunkn,nunkn*(jnode-1)+junkn,:)= squeeze(M(nunkn*(inode-1)+iunkn,nunkn*(jnode-1)+junkn,:)) ...
                                    + v(:)/dt.*geometry_variable.dvolu(:,igauss);
%                                 obj.LHS(iv,jv,:) = squeeze(obj.LHS(iv,jv,:)) + v(:).*geometry.dvolu(:,igauss);
                            end
                        end
                    end
                end
             end
            obj.M_elem = M;
        end
        
        function K = compute_K(obj,dim,nelem,geometry_test,geometry_variable,material)
             nunkn = dim.nunkn(1);
            K = zeros(nunkn*geometry_test.nnode,nunkn*geometry_variable.nnode,nelem);
           
            Cmat = material.mu;
            
            for igauss = 1 :geometry_variable.ngaus
                
                Bmat = obj.computeB(nunkn,nelem,geometry_variable.nnode,geometry_variable.cartd(:,:,:,igauss));
                
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
            obj.K_elem = K;           
            end
        
        function D = compute_D(obj,dim,nelem,geometry_test,geometry_variable)
            nunkn_u=dim.nunkn(1);
 
            
            D = zeros(nunkn_u*geometry_test.nnode,geometry_variable.nnode,nelem);

            for igauss=1:geometry_test.ngaus
                for inode_var = 1:geometry_variable.nnode
                    for inode_test = 1:geometry_test.nnode
                        for idime = 1:geometry_test.ndime
                            dof_test = inode_test*nunkn_u - nunkn_u + idime;
                            v= squeeze (geometry_test.cartd(idime,inode_test,:,igauss));
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
        
       function Fext = compute_vol_force_on_nodes(obj,geometry_variable,idx,nnode,nunkn)
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
       end
        
       function Fext = compute_vol_force_on_gauss_points(obj,geometry_variable,nnode,nunkn,f)
           Fext = zeros(nnode*nunkn,1,obj.nelem);
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
       
       function g = compute_velocity_divergence(obj)
           g = zeros(obj.geometry(2).nnode*obj.dim.nunkn(2),1,obj.nelem);
       end
       
       function variable = computeVars(obj,x_free)
            x = obj.reduced_vector_2_full_vector(x_free);
            variable.u = x(1:obj.dof.ndof(1),:);
            variable.p = x(obj.dof.ndof(1)+1:end,:);
       end
        
    end
    methods (Access=protected)
        function mat = computeMatrix(obj,dim,nelem,geometry_test,geometry_variable,material,dt,ifield,jfield)

           if ifield == 1 && jfield==1
               K = obj.compute_K(dim,nelem,geometry_test,geometry_variable,material);
               M = obj.compute_M(dim,nelem,geometry_test,geometry_variable,dt);
               mat = M + K;
           elseif ifield == 1 && jfield==2
               mat = obj.compute_D(dim,nelem,geometry_test,geometry_variable);
           elseif ifield == 2 && jfield==1
               mat = permute(obj.LHS_elem{1,2},[2,1,3]);
           else
               D = obj.LHS_elem{2,1};
               D_traspose= obj.LHS_elem{1,2};
               row = length(D(:,1,1));
               col = length (D_traspose(1,:,1));
%                D = obj.LHS_elem{1,1};
%                D_traspose= obj.LHS_elem{1,1};
%                row = length(D(:,1,1));
%                col = length (D_traspose(1,:,1));

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
        function Fext = computeSuperficialFext(obj,nunkn,nelem,nnode,bc,idx) %To be donne
%             Fext = zeros(nnode*nunkn,1,nelem);
              Fext=0;
        end
        function Fext = computeVolumetricFext(obj,dim,nelem,geometry_variable,dof)
            idx = obj.dof.in_elem{1};
            geometry_variable = geometry_variable(1);
            nnode = geometry_variable.nnode;
            nunkn= dim.nunkn(1);
%             f = zeros(nnode*nunkn,1,nelem);
           
%             obj.RHS = zeros(nnode*nunkn,1,nelem);


            if  ~isempty(dof.neumann_values)
                if ~ismatrix(dof.neumann_values) 
                    f=dof.neumann_values;
                    Fext = obj.compute_vol_force_on_gauss_points(geometry_variable,nnode,nunkn,f);
                else
                    Fext = obj.compute_vol_force_on_nodes(geometry_variable,idx,nnode,nunkn);
                end
            else
                 Fext = zeros(nnode*nunkn,1,nelem);
            end

        end
    end
    
end