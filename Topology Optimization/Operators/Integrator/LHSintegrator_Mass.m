classdef LHSintegrator_Mass < LHSintegrator

    properties (Access = public)
        elemMass
    end
    
    properties (Access = private)
        quadType
    end

    
    methods (Access = public)
        
        function obj = LHSintegrator_Mass(cParams)
            obj.init(cParams);
            obj.quadType = cParams.quadType;
            obj.createQuadrature();
            obj.createInterpolation();
            if isfield(cParams, 'interpolation')
                obj.interpolation = cParams.interpolation;
                obj.quadrature    = cParams.quadrature;
            end
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);
        end

        function lhs = computeElemental(obj) % Temporary for Stokes
            lhs = obj.computeElementalLHS();
        end
        
    end
    
    methods (Access = protected)
        
        function lhs = computeElementalLHS(obj)
            shapes = obj.interpolation.shape;
            quad   = obj.quadrature;
            dvolu  = obj.mesh.computeDvolume(quad);
            ngaus  = obj.quadrature.ngaus;
            nelem  = obj.mesh.nelem;
%             nnode  = obj.mesh.nnodeElem;
            ndimf  = obj.dim.ndimf;
            nnode  = obj.interpolation.nnode;

            % One dimension
%             lhs = zeros(nnode,nnode,nelem);
%             for igaus = 1:ngaus
%                 dv(1,1,:) = dvolu(igaus,:);
%                 Ni = shapes(:,igaus);
%                 Nj = shapes(:,igaus);
%                 NiNj = Ni*Nj';
%                 Aij = bsxfun(@times,NiNj,dv);
%                 lhs = lhs + Aij;
%             end

            % N dimensions, pending optimization
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
                                v = squeeze(Ni.*Nj);
                                M(idof, jdof, :)= squeeze(M(idof,jdof,:)) ...
                                    + v(:).*dvol;
                            end
                        end
                    end
                end
            end
            lhs = M;
            obj.elemMass = lhs;  
        end

       function createQuadrature(obj)
           quad = Quadrature.set(obj.mesh.type);
           quad.computeQuadrature(obj.quadType);
           obj.quadrature = quad;
       end
        
    end
    
end