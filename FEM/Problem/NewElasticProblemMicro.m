classdef NewElasticProblemMicro < NewElasticProblem %NewFEM

    properties (Access = private)
        variables2print
    end

    methods (Access = public)

        function setMatProps(obj,s)
           obj.material.compute(s);
        end

        function mesh = getMesh(obj)
            mesh  = obj.mesh;
        end
        
        function interp = getInterpolation(obj)
            interp  = obj.interp{1};
        end

        function quad = getQuadrature(obj)
            quad  = obj.quadrature;
        end
        
        function v = computeGeometricalVolume(obj)
            v = 1;%sum(sum(obj.geometry.dvolu));
        end

        function Ch = computeChomog(obj)
            nstre = obj.dim.nstre;
            ngaus = obj.dim.ngaus;
            nelem = obj.dim.nelem;
            ndof  = obj.dim.ndof;
            basis = diag(ones(nstre,1));
            tStrn  = zeros(nstre,ngaus,nstre,nelem);
            tStrss = zeros(nstre,ngaus,nstre,nelem);
            tDisp  = zeros(nstre,ndof);
            Ch = zeros(nstre,nstre);
            for istre=1:nstre
                obj.vstrain = basis(istre,:);
                obj.solve();
                vars = obj.computeStressStrainAndCh();
                Ch(:,istre)         = vars.stress_homog;
                tStrn(istre,:,:,:)  = vars.strain;
                tStrss(istre,:,:,:) = vars.stress;
                tDisp(istre,:)      = vars.d_u;
            end
            obj.variables.Chomog  = Ch;
            obj.variables.tstrain = tStrn;
            obj.variables.tstress = tStrss;
            obj.variables.tdisp   = tDisp;
        end

    end

    methods (Access = private)

        function vars = computeStressStrainAndCh(obj)
            vars = obj.variables;
            vars.stress_fluct = vars.stress;
            vars.strain_fluct = vars.strain;
            Cmat  = obj.material.C;
            ngaus = obj.dim.ngaus;
            nstre = obj.dim.nstre;
            nelem = obj.dim.nelem;
            dV = obj.getDvolume()';
            
            vars.stress = zeros(ngaus,nstre,nelem);
            vars.strain = zeros(ngaus,nstre,nelem);
            vars.stress_homog = zeros(nstre,1);
            vol_dom = 1;%sum(sum(obj.geometry.dvolu));
            
            for igaus = 1:ngaus
                vars.strain(igaus,1:nstre,:) = obj.vstrain.*ones(1,nstre,nelem) + vars.strain_fluct(igaus,1:nstre,:);
                for istre = 1:nstre
                    for jstre = 1:nstre
                        Cij = squeeze(Cmat(istre,jstre,:,igaus));
                        C  = squeeze(Cij);
                        strs = squeeze(vars.stress(igaus,istre,:));
                        strn = squeeze(vars.strain(igaus,jstre,:));
                        vars.stress(igaus,istre,:) = strs + 1/vol_dom* C.* strn;
                    end
                end
                % Contribution to Chomogenized
                for istre = 1:nstre
                    strs = squeeze(vars.stress(igaus,istre,:));
                    vars.stress_homog(istre) = vars.stress_homog(istre) +  1/vol_dom *(strs)'*dV(:,igaus);
                end
                
            end
            obj.variables = vars;
        end

    end

end