classdef LHSintegrator_MassModal < LHSintegrator
    
    properties (Access = public)
        elemMass
    end
    
    properties (Access = private)
        quadType
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = LHSintegrator_MassModal(cParams)
            obj.init(cParams);
            obj.quadType = cParams.quadType;
            obj.createQuadrature();
            obj.createInterpolation();
            if isfield(cParams, 'interpolation')
                obj.interpolation = cParams.interpolation;
                obj.quadrature    = cParams.quadrature;
            end
        end

        function LHS = compute(obj,x)
            lhs = obj.computeElementalLHS(x);
            LHS = obj.assembleMatrix(lhs);
        end

%         function lhs = computeElemental(obj)
%             lhs = obj.computeElementalLHS();
%         end
        
    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj,x)
            rho    = x{1};
            shapes = obj.interpolation.shape;
            quad   = obj.quadrature;
            dvolu  = obj.mesh.computeDvolume(quad);
            ngaus  = obj.quadrature.ngaus;
            nelem  = obj.mesh.nelem;
            ndimf  = obj.dim.ndimf;
            nnode  = obj.interpolation.nnode;
            M = zeros(nnode*ndimf,nnode*ndimf,nelem);
            dvolu = dvolu';
            for igauss = 1 :ngaus
                for inode= 1:nnode
                    for jnode= 1:nnode
                        for iunkn= 1:ndimf
                            for junkn= 1:ndimf
                                idof = ndimf*(inode-1)+iunkn;
                                jdof = ndimf*(jnode-1)+junkn;
                                dvol = dvolu(:,igauss);
                                Ni = shapes(inode,igauss,:);
                                Nj = shapes(jnode,igauss,:);
                                v = squeeze(Ni.*Nj).*rho;
                                M(idof, jdof, :)= squeeze(M(idof,jdof,:)) ...
                                    + v(:,igauss).*dvol;
                            end
                        end
                    end
                end
            end
            lhs = M;
            obj.elemMass = lhs;
            lhs = obj.addNonStructuralMass(lhs);
        end

        function lhs = addNonStructuralMass(obj,lhs)
            %% CantileverBeam_Triangular
            % Center of the right side
%                         lhs(:,:,2) = lhs(:,:,2) + 10*eye(6,6);
%                         lhs(:,:,7) = lhs(:,:,7) + 10*eye(6,6);
            % Corners of the right side
            %             lhs(:,:,214) = lhs(:,:,214) + 20*lhs(:,:,214);
            %             lhs(:,:,8) = lhs(:,:,8) + 10*lhs(:,:,8);
            %% CantileverArnau2
            lhs(:,:,1502) = lhs(:,:,1502) + 1*eye(6,6);
 %           lhs(:,:,3022) = lhs(:,:,3022) + 10*eye(6,6);
 %           lhs(:,:,1500) = lhs(:,:,1500) + 10*eye(6,6);
 %           lhs(:,:,1503) = lhs(:,:,1503) + 10*eye(6,6);
%             lhs(:,:,1461) = lhs(:,:,1461) + 10*eye(6,6);
%             lhs(:,:,1462) = lhs(:,:,1462) + 10*eye(6,6);
%             lhs(:,:,1463) = lhs(:,:,1463) + 10*eye(6,6);
%             lhs(:,:,1464) = lhs(:,:,1464) + 10*eye(6,6);
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature(obj.quadType);
            obj.quadrature = quad;
        end

    end

%     methods (Access = protected)
% 
%         function M = assembleMatrix(obj,Melem,x)
%             M = obj.assembleMatrixViaIndices(Melem,x);
%         end
% 
%     end
% 
%     methods (Access = private)
% 
%         function M = assembleMatrixViaIndices(obj,Melem,xReg)
%              MeRe = obj.restricMelem(Melem,xReg);
%              connec    = obj.globalConnec;
%              nnodes  = obj.dim.nnodes;
%              ndimf   = obj.dim.ndimf;
%              ndofs   = ndimf*nnodes;
%              nelem   = size(connec, 1);
%              dofConnec = obj.computeDofConnectivity()';
%              ndofEl  = size(dofConnec,2);
%              res = zeros(ndofEl^2 * nelem, 3);
%              strt = 1;
%              fnsh = nelem;
%              ndofEl1 = size(MeRe,1);
%              ndofEl2 = size(MeRe,2);
%              for i = 1:ndofEl1
%                  dofsI = dofConnec(:,i);
%                  for j = 1:ndofEl2
%                      dofsJ = dofConnec(:,j);
%                      a = squeeze(MeRe(i,j,:));
%                      matRes = [dofsI, dofsJ, a];
%                      res(strt:fnsh,:) = matRes;
%                      strt = strt + nelem;
%                      fnsh = fnsh + nelem;
%                  end
%              end
%              M = sparse(res(:,1), res(:,2), res(:,3), ndofs, ndofs);
%         end
% 
%         function MeRe = restricMelem(obj,Melem,xReg)
%             xReg = xReg{1};
%             nelem = size(xReg,1);
%             MeRe = zeros(size(Melem));
%             c1 = 6*10^5;
%             c2 = -5*10^6;
%             for iElem = 1: nelem
%                 if xReg(iElem,1) >= 0.1
%                     MeRe(:,:,iElem) = Melem(:,:,iElem); % xReg(iElem)*
%                 else
%                     MeRe(:,:,iElem) = (c1*(xReg(iElem))^4+c2*(xReg(iElem))^5)*Melem(:,:,iElem);
%                 end
%             end
%         end
% 
%         function dofConnec = computeDofConnectivity(obj)
%             connec  = obj.globalConnec;
%             ndimf   = obj.dim.ndimf;
%             nnodeEl = size(connec, 2); % obj.dim.nnodeElem
%             ndofsEl = nnodeEl * ndimf; % obj.dim.ndofsElem;
%             dofsElem  = zeros(ndofsEl,size(connec,1));
%             for inode = 1:nnodeEl
%                 for iunkn = 1:ndimf
%                     idofElem   = obj.nod2dof(inode,iunkn);
%                     globalNode = connec(:,inode);
%                     idofGlobal = obj.nod2dof(globalNode,iunkn);
%                     dofsElem(idofElem,:) = idofGlobal;
%                 end
%             end
%             dofConnec = dofsElem;
%         end
% 
%         function idof = nod2dof(obj, inode, iunkn)
%             ndimf = obj.dim.ndimf;
%             idof(:,1)= ndimf*(inode - 1) + iunkn;
%         end

end
