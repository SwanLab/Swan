classdef ElasticProblemMicro < ElasticProblem

    properties (Access = public)
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
                v2p = obj.assignVarsToPrint(istre);
            end
            obj.variables2print   = v2p;
            obj.variables.Chomog  = Ch;
            obj.variables.tstrain = tStrn;
            obj.variables.tstress = tStrss;
            obj.variables.tdisp   = tDisp;
        end

    end

    methods (Access = private)

        function vars = computeStressStrainAndCh(obj)
            vStrn = obj.vstrain;
            vars  = obj.variables;
            Cmat  = obj.material.C;
            ngaus = obj.dim.ngaus;
            nstre = obj.dim.nstre;
            nelem = obj.dim.nelem;
            strFluct = vars.strain;
            dV = obj.getDvolume()';
            
            vars.stress_fluct = vars.stress;
            vars.strain_fluct = vars.strain;
            
            vars.stress = zeros(ngaus,nstre,nelem);
            vars.strain = zeros(ngaus,nstre,nelem);
            vars.stress_homog = zeros(nstre,1);
            
            for igaus = 1:ngaus
                vars.strain(igaus,1:nstre,:) = vStrn.*ones(1,nstre,nelem) + strFluct(igaus,1:nstre,:);
                for istre = 1:nstre
                    for jstre = 1:nstre
                        Cij  = squeeze(Cmat(istre,jstre,:,igaus));
                        C    = squeeze(Cij);
                        strs = squeeze(vars.stress(igaus,istre,:));
                        strn = squeeze(vars.strain(igaus,jstre,:));
                        vars.stress(igaus,istre,:) = strs + C.* strn;
                    end
                    strs = squeeze(vars.stress(igaus,istre,:));
                    vars.stress_homog(istre) = vars.stress_homog(istre) + (strs)'*dV(:,igaus);
                end
            end
            obj.variables = vars;
        end

        function v2p = assignVarsToPrint(obj, istre)
            vars = obj.variables;
            v2p{istre}.d_u          = vars.d_u;
            v2p{istre}.fext         = vars.fext;
            v2p{istre}.stress       = vars.stress;
            v2p{istre}.strain       = vars.strain;
            v2p{istre}.stress_fluct = vars.stress_fluct;
            v2p{istre}.strain_fluct = vars.strain_fluct;
        end
    end

end