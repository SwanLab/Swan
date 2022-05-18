classdef Element_Stokes < Element
    %Element_Stokes Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        LHS_elem
        LHS
        RHS
    end

    properties(Access = private)
        D_elem
        M_elem
        K_elem
        D
        dt
        mesh
        velocityField
        pressureField
        forcesFormula
    end
    
    methods
        function obj = Element_Stokes(geometry,mesh,material,dof,problemData,interp, vField, pField,forcesFormula)
            obj.initElement(geometry,mesh,material,dof,problemData.scale,interp);
%             obj.dim = dim;
            obj.mesh = mesh;
            %obj.nstre=0;
            obj.nfields=2;
            obj.velocityField = vField;
            obj.pressureField = pField;
            obj.forcesFormula = forcesFormula;
        end
        
        function [r,dr] = computeResidual(obj,x,dr,x_n)
%             K = compute_LHS(obj);
            freeV = obj.velocityField.boundaryConditions.free;
            lenFreeV = length(freeV);
            if (nargin ==3)
                % Steady
                Mred_x_n = zeros(lenFreeV,1);
            else
                % Transient
            
                s.dim          = obj.velocityField.dim;
                s.globalConnec = obj.velocityField.connec;
                s.nnodeEl      = obj.velocityField.dim.nnodeElem;
                assembler = Assembler(s);
                M = assembler.assemble(obj.M_elem);
                M = obj.symGradient(M);

%                 M = obj.AssembleMatrix(obj.M_elem,1,1);
%                 M = obj.symGradient(M);
                Mred = M(freeV,freeV);
                Mred_x_n = Mred*x_n;
            end
            
            Fext = compute_RHS(obj);
            
            
            R = obj.compute_imposed_displacement_force(obj.LHS);
            Fext = Fext + R ;
            
            
            Fext_red = obj.bcApplier.fullToReducedVector(Fext);
            Fext_red(1:lenFreeV,1) = Fext_red(1:lenFreeV,1) + Mred_x_n;
            
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
            obj.dt = dt;
            AA = obj.computeVelocityLaplacian();
            D = obj.computeDmatrix();
            BB = obj.computePressureLHSMatrix();
            LHS = [AA, D; D',BB];
        end
        
        function RHS = compute_RHS(obj)
            % Inefficient. It is always the same.
            Fext = obj.computeVolumetricFext();
            g = obj.compute_velocity_divergence;
            RHS_elem{1,1} = Fext;
            RHS_elem{2,1} = g;
            RHS = AssembleVector(obj,RHS_elem);
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
        
        function Fext = compute_vol_force_on_gauss_points(obj)
            geometry = obj.velocityField.geometry;
            shapesV  = obj.velocityField.interpolation.shape;
            dvol = geometry.dvolu;
            ngaus = size(dvol,2);
            nnode = obj.velocityField.dim.nnodeElem;
            nunkn = obj.velocityField.dim.ndimf;
            Fext = zeros(nnode*nunkn,1,obj.nelem);
            
            f = obj.calculateForcesFromExpression();
                
            for igaus=1:ngaus
                for inode=1:nnode
                    for iunkn=1:nunkn
                        elemental_dof = inode*nunkn-nunkn+iunkn; %% dof per guardar el valor de la integral
                        shape = shapesV(inode,igaus);
                        fvalue = f(iunkn,igaus,:);
                        v= squeeze(shape.*fvalue);
                        Fext(elemental_dof,1,:)= squeeze(Fext(elemental_dof,1,:)) + v(:).*dvol(:,igaus);
                        
                    end
                end
            end
        end

        function f = calculateForcesFromExpression(obj)
            ngaus  = size(obj.velocityField.interpolation.shape,2);
            xGauss = obj.velocityField.xGauss;
            nelem = obj.mesh.nelem;
            for ielem = 1:nelem
                ind=1;
                for igaus = 1:ngaus
                    xG = xGauss(:,igaus,ielem);
                    pos_node = num2cell(xG);
                    fCell = obj.forcesFormula(pos_node{:});
                    fMat = cell2mat(fCell);
                    F(:,igaus,ielem) = fMat;
                    ind=ind+length(fMat);
                end
            end
            f = F;
        end
        
        function g = compute_velocity_divergence(obj)
            nunkn = obj.velocityField.dim.ndimf;
            nnode = obj.pressureField.dim.nnodeElem;
            g = zeros(nnode*nunkn,1,obj.nelem);
        end
        
        function variable = computeVars(obj,x_free)
            x = obj.bcApplier.reducedToFullVector(x_free);
            ndofsV = obj.velocityField.dim.ndofs;
            variable.u = x(1:ndofsV,:);
            variable.p = x(ndofsV+1:end,:);
        end
    end
    
    methods (Access = protected)
        
%         function Fext = computePuntualRHS(obj,nunkn,nelem,nnode,bc,idx)
%             Fext = zeros(nnode*nunkn,1,nelem);
%             for i = 1:length(bc.iN)
%                 for j = 1:nelem
%                     ind = find(idx(:,j) == bc.iN(i));
%                     if ~isempty(ind)
%                         Fext(ind,:,j) = bc.neunodes(i,3);
%                     end
%                     % clear ind
%                     ind = [];
%                 end
%             end
%         end
        
%         function Fext = computeSuperficialFext(obj,nunkn,nelem,nnode,bc,idx) %To be donne
%             % Fext = zeros(nnode*nunkn,1,nelem);
%             Fext = 0;
%         end
        
        function Fext = computeVolumetricFext(obj)
            geometry = obj.velocityField.geometry;
            nnode = obj.velocityField.dim.nnodeElem;
            nunkn = obj.velocityField.dim.ndimf;
            nelem = obj.nelem;
            Fext = obj.compute_vol_force_on_gauss_points();
%             f = zeros(nnode*nunkn,1,nelem);
%             obj.RHS = zeros(nnode*nunkn,1,nelem);
            
%             if  ~isempty(dof.neumann_values)
%                 if ~ismatrix(dof.neumann_values)
%                     Fext = obj.compute_vol_force_on_gauss_points();
%                 else
% %                     idx = obj.dof.in_elem{1};
% %                     Fext = obj.compute_vol_force_on_nodes(geometry,idx,nnode,nunkn);
%                 end
%             else
% %                 Fext = zeros(nnode*nunkn,1,nelem);
%             end
        end
    end

    methods (Access = private)

        function A = symGradient(obj, B)
            A = 1/2 * (B+B');
        end

        function Aassembled = computeVelocityLaplacian(obj)
            obj.computeMelem();
            obj.computeKelem();
            A = obj.K_elem + obj.M_elem;
            
            s.dim          = obj.velocityField.dim;
            s.globalConnec = obj.velocityField.connec;
            s.nnodeEl      = obj.velocityField.dim.nnodeElem;
            assembler = Assembler(s);
            lhs = assembler.assemble(A);
            Aassembled = obj.symGradient(lhs);

%             Aass = obj.AssembleMatrix(A ,1, 1);
%             Aassembled = obj.symGradient(Aass);
        end

        function Dassembled = computeDmatrix(obj)
            obj.computeDelem();
%             Dassembled = obj.AssembleMatrix(obj.D_elem, 1, 2);
%             obj.D = Dassembled;
            s.dim = [];
            s.nnodeEl = [];
            s.globalConnec = [];
            assembler = Assembler(s);
            Dassembled = assembler.assembleFields(obj.D_elem, ...
                                    obj.velocityField, obj.pressureField);
            obj.D = Dassembled;
        end

        function BB = computePressureLHSMatrix(obj)
            sz = size(obj.D, 2);
            BB = sparse(sz,sz);
        end
        
        function computeMelem(obj)
            vel = obj.velocityField;
            nunkn = vel.dim.ndimf;
            nnode = vel.dim.nnodeElem;
            ndofs = vel.dim.ndofsElem;
            nelem = obj.nelem;
            dtime = obj.dt;
            shpeV = vel.interpolation.shape;
            dvolV = vel.geometry.dvolu;
            ngaus = size(dvolV,2);
            M = zeros(ndofs, ndofs, nelem);
            
            for igauss = 1 :ngaus
                for inode= 1:nnode
                    for jnode= 1:nnode
                        for iunkn= 1:nunkn
                            for junkn= 1:nunkn
                                idof = nunkn*(inode-1)+iunkn;
                                jdof = nunkn*(jnode-1)+junkn;
                                dvol = dvolV(:,igauss);
                                Ni = shpeV(inode,igauss,:);
                                Nj = shpeV(jnode,igauss,:);
                                v = squeeze(Ni.*Nj);
                                M(idof, jdof, :)= squeeze(M(idof,jdof,:)) ...
                                    + v(:)/dtime.*dvol;
                            end
                        end
                    end
                end
            end
            
            obj.M_elem = M;
%             obj.computeMassMatrix();
        end
        
        function computeKelem(obj)
            vel = obj.velocityField;
            ndofs = vel.dim.ndofsElem;
            nelem = obj.nelem;

            material = obj.material;
            Cmat = material.mu;
            
            geom = vel.geometry;
            shape = geom.dNdx;
            ngaus = size(shape,4);
            dvolu = geom.dvolu;
            lhs = zeros(ndofs, ndofs, nelem);
            for igaus = 1:ngaus
                dNdx = shape(:,:,:,igaus);
                dV(1,1,:) = dvolu(:,igaus);
                Bmat = obj.computeB(dNdx);
                Bt   = permute(Bmat,[2 1 3]);
%                 BtC  = pagemtimes(Bt,Cmat);
%                 BtCB = pagemtimes(BtC, Bmat);
                BtB = pagemtimes(Bt,Bmat);
                lhs = lhs + bsxfun(@times, BtB, dV);
            end
            obj.K_elem = lhs;
        end
        
        function computeDelem(obj)
            nelem = obj.nelem;
            vel = obj.velocityField;
            prs = obj.pressureField;
            nunknV = vel.dim.ndimf;
            nnodeV = vel.dim.nnodeElem;
            nnodeP = prs.dim.nnodeElem;
            
            dNdxV = vel.geometry.dNdx;
            dvolV = vel.geometry.dvolu;
            shpeP = prs.interpolation.shape; %nope, should be Quadratic
            ngaus = size(dNdxV,4);

            D = zeros(nunknV*nnodeV,nnodeP,nelem);
            for igauss=1:ngaus
                for inode_var = 1:nnodeP
                    for inode_test = 1:nnodeV
                        for idime = 1:vel.interpolation.ndime
                            dof_test = inode_test*nunknV - nunknV + idime;
                            v = squeeze(dNdxV(idime,inode_test,:,igauss));
                            D(dof_test,inode_var,:)= squeeze(D(dof_test,inode_var,:)) - v(:).*shpeP(inode_var,igauss)...
                                .*dvolV(:,igauss);
                        end
                    end
                end
            end
            obj.D_elem = D;
        end
        
        function B = computeB(obj,dNdx)
            vel = obj.velocityField;
            nunkn = vel.dim.ndimf;
            nnode = vel.dim.nnodeElem;
            nelem = obj.nelem;
            B = zeros(2,nnode*nunkn,nelem); %check the 2
            for i = 1:nnode
                j = nunkn*(i-1)+1;
                B(1,j,:)  = dNdx(1,i,:);
                B(2,j+1,:)= dNdx(1,i,:);
                B(3,j,:)  = dNdx(2,i,:);
                B(4,j+1,:)= dNdx(2,i,:);
            end
        end
    end

end