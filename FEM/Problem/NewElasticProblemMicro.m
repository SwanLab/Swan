classdef NewElasticProblemMicro < NewElasticProblem %NewFEM

    properties (Access = private)
        Chomog
        variables2print
        tstrain
        tstress
    end

    methods (Access = public)

        function setMatProps(obj,s)
           obj.material.compute(s);
        end

        function solveMicro(obj)
            obj.computeStressStrainAndCh();
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

        function [Ch,tstrain,tstress] = computeChomog(obj)
            nstre = obj.dim.nstre;
            ngaus = obj.dim.ngaus;
            nelem = obj.dim.nelem;
            ndof  = obj.dim.ndof;
            basis = diag(ones(nstre,1));
            tstrain  = zeros(nstre,ngaus,nstre,nelem);
            tstrainF = zeros(nstre,ngaus,nstre,nelem);
            tstress  = zeros(nstre,ngaus,nstre,nelem);
            tdisp    = zeros(nstre,ndof);
            
            Ch2 = zeros(nstre,nstre);
            Ch = zeros(nstre,nstre);
            
            var2print = cell(nstre,1);
            for istre=1:nstre
                obj.vstrain = basis(istre,:);
                obj.solve();
                obj.computeStressStrainAndCh();
                Ch(:,istre)           = obj.variables.stress_homog;
                tstrain(istre,:,:,:)  = obj.variables.strain;
                tstrainF(istre,:,:,:) = obj.variables.strain_fluct;
                tstress(istre,:,:,:)  = obj.variables.stress;
                tdisp(istre,:)        = obj.variables.d_u;
                var2print{istre} = obj.assignVarsToPrint();
            end
            obj.variables2print = var2print;

            obj.variables.Chomog  = Ch;
            obj.variables.tstrain = tstrain;
            obj.variables.tstress = tstress;
            obj.variables.tdisp   = tdisp;
            
            dV = obj.getDvolume()';
            for istre = 1:nstre
                for jstre = 1:nstre
                    s = squeezeParticular(tstress(istre,:,:,:),1);
                    e = squeezeParticular(tstrain(jstre,:,:,:),1);
                    ener = (s.*e);
                    en = sum(ener,2);
                    en = squeezeParticular(en,2);
                    Ch2(istre,jstre) = en(:)'*dV(:);
                end
            end
            
            
            Cmat = obj.material.C;
            
            Ch3 = zeros(nstre,nstre);
            ngaus = size(tstrain,2);
            nelem = size(tstrain,4);
            for istre = 1:nstre
                for jstre = 1:nstre
                   ei(1:ngaus,:,:) = squeeze(tstrain(istre,:,:,:));
                   ej(1:ngaus,:,:) = squeeze(tstrain(jstre,:,:,:));
                   c = zeros(ngaus,nelem);
                    for kstre = 1:nstre
                        for lstre = 1:nstre
                         eiV(1:ngaus,:) = squeeze(ei(:,kstre,:));
                         ejV(1:ngaus,:) = squeeze(ej(:,lstre,:));
                         Cm = squeezeParticular(Cmat(kstre,lstre,:,:),1);
                         Cm = squeezeParticular(Cm,1)';
                         %Cm  = repmat(Cmm,ngaus,1);
                         c = c + Cm.*eiV.*ejV;
                        end
                    end
                    cC = c.*dV';
                    Ch3(istre,jstre) = sum(cC(:));
                end
            end
            
            obj.Chomog = Ch;
            obj.tstrain = tstrain;
            obj.tstress = tstress;
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

        function microVars = assignMicroVars(obj, istre)
        end

        function v2p = assignVarsToPrint(obj)
            vars = obj.variables;
            v2p.d_u          = vars.d_u;
            v2p.fext         = vars.fext;
            v2p.stress       = vars.stress;
            v2p.strain       = vars.strain;
            v2p.stress_fluct = vars.strain_fluct;
            v2p.strain_fluct = vars.strain_fluct;
        end

    end

end