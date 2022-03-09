classdef NewElasticProblemMicro < NewElasticProblem %NewFEM

    properties (Access = private)
        Chomog
        variables2print
        tstrain
        tstress
        vstrain
    end

    methods (Access = public)

        function solveMicro(obj)
            obj.computeStressStrainAndCh();
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
                strain = basis(istre,:);
                obj.vstrain = strain;
                obj.computeStressStrainAndCh();
                obj.solve();
                Ch(:,istre) = obj.variables.stress_homog;
                tstrain(istre,:,:,:) = obj.variables.strain;
                tstrainF(istre,:,:,:) = obj.variables.strain_fluct;
                tstress(istre,:,:,:) = obj.variables.stress;
                
                tdisp(istre,:) = obj.variables.d_u;
                var2print{istre}.stress = obj.variables.stress;
                var2print{istre}.strain = obj.variables.strain;
                var2print{istre}.stress_fluct = obj.variables.strain_fluct;
                var2print{istre}.strain_fluct = obj.variables.strain_fluct;
                var2print{istre}.d_u = obj.variables.d_u;
                var2print{istre}.fext = obj.variables.fext;
            end
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
            

            
            obj.variables2print = var2print;
            
            obj.Chomog = Ch;
            obj.tstrain = tstrain;
            obj.tstress = tstress;
        end

    end

    methods (Access = private)

        function vars = computeStressStrainAndCh(obj)
            vars = obj.variables;
            vars.stress_fluct = permute(vars.stress,[2 3 1]);
            vars.strain_fluct = permute(vars.strain,[2 3 1]);
            Cmat = obj.material.C;
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
                        C  = squeeze(squeeze(Cmat(istre,jstre,:,igaus)));
                        vars.stress(igaus,istre,:) = squeeze(vars.stress(igaus,istre,:)) + 1/vol_dom*C.* squeeze(vars.strain(igaus,jstre,:));
                    end
                end
                % Contribution to Chomogenized
                for istre = 1:nstre
                    vars.stress_homog(istre) = vars.stress_homog(istre) +  1/vol_dom *(squeeze(vars.stress(igaus,istre,:)))'*dV(:,igaus);
                end
                
            end
            
            obj.variables = vars;
        end
        
    end

end