classdef LHSintegrator_StokesD < LHSintegrator

    properties (Access = private)
        pressure
        velocity
    end

    methods (Access = public)

        function obj = LHSintegrator_StokesD(cParams)
            %             obj.init(cParams);
            %             obj.createQuadrature();
            %             obj.createInterpolation();
            %             obj.createGeometry();
            obj.initStokesD(cParams);
        end

        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleStokesD(lhs);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
            vel = obj.velocity;
            prs = obj.pressure;
            nelem = obj.mesh.nelem;
            nunknV = vel.dim.ndimf;
            nnodeV = vel.dim.nnodeElem;
            nnodeP = prs.dim.nnodeElem;
            
            dNdxV = vel.geometry.dNdx;
            dvolV = vel.geometry.dvolu;
            shpeP = prs.interpolation.shape;
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
            lhs = D;
        end

    end

    methods (Access = private)

        function initStokesD(obj, cParams)
            obj.mesh     = cParams.mesh;
            obj.pressure = cParams.pressure;
            obj.velocity = cParams.velocity;
%             obj.material = cParams.material;
        end

        function LHS = assembleStokesD(obj,Delem)
            s.dim           = [];
            s.nnodeEl       = [];
            s.globalConnec  = [];
            assembler = Assembler(s);
            LHS = assembler.assembleFields(Delem,obj.velocity,obj.pressure);
        end

    end

end