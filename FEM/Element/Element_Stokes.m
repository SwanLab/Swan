classdef Element_Stokes < Element
    %Element_Stokes Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        LHS_elem
        LHS
        RHS
        M_elem
        K_elem
        interpolation_v
        interpolation_p
        dof_v
        dof_p
    end
    
    methods
        function obj = Element_Stokes(geometry,mesh,material,dof,problemData,interp)
            obj.initElement(geometry,mesh,material,dof,problemData.scale,interp);
            %obj.nstre=0;
            obj.nfields=2;
            obj.interpolation_v= interp{1};
            obj.interpolation_p= interp{2};
        end
        
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
            
            R = obj.compute_imposed_displacement_force(obj.LHS);
            Fext = Fext + R ;
            
            %             Kred = obj.fullToReducedMatrix(K);
            
            Fext_red = obj.bcApplier.fullToReducedVector(Fext);
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
            LHSred = obj.bcApplier.fullToReducedMatrix(obj.LHS);
            dr = LHSred;
        end
        
        function LHS = compute_LHS(obj,dt)
            for ifield = 1:obj.nfields
                for jfield = 1:obj.nfields
                    obj.LHS_elem{ifield,jfield} = obj.computeMatrix(obj.nelem,obj.material,dt,ifield,jfield);
                    % mat = obj.computeMatrix(obj.nelem,obj.geometry(ifield),obj.geometry(jfield),obj.material,ifield,jfield);
                    LHS{ifield,jfield} = obj.AssembleMatrix(obj.LHS_elem{ifield,jfield},ifield,jfield);
                end
            end
            % LHS= AssembleMatrix(obj,obj.LHS_elem);
            LHS = cell2mat(LHS);
        end
        
        function RHS = compute_RHS(obj)
            Fext = obj.computeVolumetricFext(obj.nelem,obj.geometry,obj.dof);
            g = obj.compute_velocity_divergence;
            RHS_elem{1,1} = Fext;
            RHS_elem{2,1} = g;
            RHS = AssembleVector(obj,RHS_elem);
        end
        
        function M = compute_M(obj,nelem,dt)
            nunkn = obj.dof.nunkn(1);
            M = zeros(nunkn*obj.interpolation_v.nnode,nunkn*obj.interpolation_v.nnode,nelem);
            
            for igauss = 1 :obj.quadrature.ngaus
                for inode= 1:obj.interpolation_v.nnode
                    for jnode= 1:obj.interpolation_v.nnode
                        for iunkn= 1:nunkn
                            for junkn= 1:nunkn
                                v = squeeze(obj.interpolation_v.shape(inode,igauss,:).*obj.interpolation_v.shape(jnode,igauss,:));
                                M(nunkn*(inode-1)+iunkn,nunkn*(jnode-1)+junkn,:)= squeeze(M(nunkn*(inode-1)+iunkn,nunkn*(jnode-1)+junkn,:)) ...
                                    + v(:)/dt.*obj.geometry(1).dvolu(:,igauss);
                                %                                 obj.LHS(iv,jv,:) = squeeze(obj.LHS(iv,jv,:)) + v(:).*geometry.dvolu(:,igauss);
                            end
                        end
                    end
                end
            end
            obj.M_elem = M;
        end
        
        function K = compute_K(obj,nelem,material)
            nunkn = obj.dof.nunkn(1);
            K = zeros(nunkn*obj.interpolation_v.nnode,nunkn*obj.interpolation_v.nnode,nelem);
            
            Cmat = material.mu;
            obj.quadrature.computeQuadrature('QUADRATIC');
            obj.geometry(1).computeGeometry(obj.quadrature,obj.interpolation_v);
            for igauss = 1 :obj.quadrature.ngaus
                Bmat = obj.computeB(nunkn,nelem,obj.interpolation_v.nnode,obj.geometry(1).cartd(:,:,:,igauss));
                %                 B_p=reshape(Bmat,[geometry.nnode*nunkn,1,nelem]);
                for iv=1:obj.interpolation_v.nnode*nunkn
                    for jv=1:obj.interpolation_v.nnode*nunkn
                        for istre=1:nunkn*obj.interpolation_v.ndime
                            for jstre=1:nunkn*obj.interpolation_v.ndime
                                v = squeeze(Bmat(istre,iv,:).*Cmat(istre,jstre,:).*Bmat(jstre,jv,:));
                                K(iv,jv,:)= squeeze(K(iv,jv,:)) + v(:).*obj.geometry(1).dvolu(:,igauss);
                                %                                 obj.LHS(iv,jv,:) = squeeze(obj.LHS(iv,jv,:)) + v(:).*geometry.dvolu(:,igauss);
                            end
                        end
                    end
                end
                
            end
            obj.K_elem = K;
        end
        
        function D = compute_D(obj,nelem)
            nunkn_u=obj.dof.nunkn(1);
            
            D = zeros(nunkn_u*obj.interpolation_v.nnode,obj.interpolation_p.nnode,nelem);
            obj.quadrature.computeQuadrature('QUADRATIC');
            obj.geometry(2).computeGeometry(obj.quadrature,obj.interpolation_p);
            for igauss=1:obj.quadrature.ngaus
                for inode_var = 1:obj.interpolation_p.nnode
                    for inode_test = 1:obj.interpolation_v.nnode
                        for idime = 1:obj.interpolation_v.ndime
                            dof_test = inode_test*nunkn_u - nunkn_u + idime;
                            v= squeeze (obj.geometry(1).cartd(idime,inode_test,:,igauss));
                            D(dof_test,inode_var,:)= squeeze(D(dof_test,inode_var,:)) - v(:).*obj.interpolation_p.shape(inode_var,igauss)...
                                .*obj.geometry(1).dvolu(:,igauss);
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
        
        function Fext = compute_vol_force_on_nodes(obj,geometry,idx,nnode,nunkn)
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
            
            
            %              for igaus=1:geometry.ngaus
            %                 for inode=1:nnode
            %                     for jnode=1:nnode
            %                         for iunkn=1:nunkn
            %                             elemental_dof = jnode*nunkn-nunkn+iunkn; %% dof per guardar el valor de la integral
            %
            %                                 v= squeeze(geometry.shape(inode,igaus).*geometry.shape(jnode,igaus).*f(elemental_dof,1,:));
            %                                 Fext(elemental_dof,1,:)= squeeze(Fext(elemental_dof,1,:)) + v(:).*geometry.dvolu(:,igaus);
            %
            %                         end
            %                     end
            %                 end
            %             end
        end
        
        function Fext = compute_vol_force_on_gauss_points(obj,geometry,nnode,nunkn,f)
            Fext = zeros(nnode*nunkn,1,obj.nelem);
            for igaus=1:obj.quadrature.ngaus
                for inode=1:nnode
                    for iunkn=1:nunkn
                        elemental_dof = inode*nunkn-nunkn+iunkn; %% dof per guardar el valor de la integral
                        shape = obj.interpolation_v.shape(inode,igaus);
                        fvalue = f(iunkn,igaus,:);
                        v= squeeze(shape.*fvalue);
                        Fext(elemental_dof,1,:)= squeeze(Fext(elemental_dof,1,:)) + v(:).*geometry.dvolu(:,igaus);
                        
                    end
                end
            end
        end
        
        function g = compute_velocity_divergence(obj)
            g = zeros(obj.interp{2}.nnode*obj.dof.nunkn(2),1,obj.nelem);
        end
        
        function variable = computeVars(obj,x_free)
            x = obj.bcApplier.reducedToFullVector(x_free);
            variable.u = x(1:obj.dof.ndof(1),:);
            variable.p = x(obj.dof.ndof(1)+1:end,:);
        end    
    end
    
    methods (Access = protected)
        function mat = computeMatrix(obj,nelem,material,dt,ifield,jfield)
            if ifield == 1 && jfield==1
                K = obj.compute_K(nelem,material);
                M = obj.compute_M(nelem,dt);
                mat = M + K;
            elseif ifield == 1 && jfield==2
                mat = obj.compute_D(nelem);
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
            % Fext = zeros(nnode*nunkn,1,nelem);
            Fext = 0;
        end
        
        function Fext = computeVolumetricFext(obj,nelem,geometry,dof)
            idx = obj.dof.in_elem{1};
            geometry = geometry(1);
            nnode = obj.interpolation_v.nnode;
            nunkn= obj.dof.nunkn(1);
            %             f = zeros(nnode*nunkn,1,nelem);
            
            %             obj.RHS = zeros(nnode*nunkn,1,nelem);
            
            if  ~isempty(dof.neumann_values)
                if ~ismatrix(dof.neumann_values)
                    f=dof.neumann_values;
                    Fext = obj.compute_vol_force_on_gauss_points(geometry,nnode,nunkn,f);
                else
                    Fext = obj.compute_vol_force_on_nodes(geometry,idx,nnode,nunkn);
                end
            else
                Fext = zeros(nnode*nunkn,1,nelem);
            end
        end
    end
end