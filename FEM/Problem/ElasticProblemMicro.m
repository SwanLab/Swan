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
            interp  = obj.mesh.interpolation;
            interp.computeShapeDeriv(obj.quadrature.posgp);
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
            v2p = cell(1,nstre);
            for istre=1:nstre
                obj.vstrain = basis(istre,:);
                obj.solve();
                vars = obj.computeStressStrainAndCh();
                Ch(:,istre)         = vars.stress_homog;
                tStrn(istre,:,:,:)  = vars.strain;
                tStrss(istre,:,:,:) = vars.stress;
                tDisp(istre,:)      = vars.d_u;
                obj.assignVarsToPrint(istre);
            end
            obj.variables.Chomog  = Ch;
            obj.variables.tstrain = tStrn;
            obj.variables.tstress = tStrss;
            obj.variables.tdisp   = tDisp;
        end

%         function print(obj,filename)
%             s.quad = obj.quadrature;
%             s.mesh = obj.mesh;
%             s.iter = 0;
%             s.variables = obj.variables2print;
%             s.ptype     = obj.problemData.ptype;
%             s.ndim      = obj.dim.ndim;
%             s.pdim      = obj.problemData.pdim;
%             s.type      = 'HomogenizedTensor';
%             fPrinter = FemPrinter(s);
%             fPrinter.print(filename);
%          end        

    end

    methods (Access = private)

        function vars = computeStressStrainAndCh(obj)
            vStrn = obj.vstrain;
            vars  = obj.variables;
            Cmat  = obj.material.C;
            ngaus = obj.dim.ngaus;
            nstre = obj.dim.nstre;
            nelem = obj.dim.nelem;
            dV = obj.mesh.computeDvolume(obj.quadrature)';
            strainFluct = vars.strain;
            stressFluct = vars.stress;
            
            stress = zeros(ngaus,nstre,nelem);
            strain = zeros(ngaus,nstre,nelem);
            stressHomog = zeros(nstre,1);
            
            for igaus = 1:ngaus
                strain(igaus,1:nstre,:) = vStrn.*ones(1,nstre,nelem) + strainFluct(igaus,1:nstre,:);
                for istre = 1:nstre
                    for jstre = 1:nstre
                        Cij  = squeeze(Cmat(istre,jstre,:,igaus));
                        C    = squeeze(Cij);
                        strs = squeeze(stress(igaus,istre,:));
                        strn = squeeze(strain(igaus,jstre,:));
                        stress(igaus,istre,:) = strs + C.* strn;
                    end
                    strs = squeeze(stress(igaus,istre,:));
                    stressHomog(istre) = stressHomog(istre) + (strs)'*dV(:,igaus);
                end
            end

            vars.stress_fluct = stressFluct;
            vars.strain_fluct = strainFluct;

            vars.stress = stress;
            vars.strain = strain;
            vars.stress_homog = stressHomog;
            obj.variables = vars;
        end

        function assignVarsToPrint(obj, istre)
            vars = obj.variables;
            ndimField = obj.dim.ndimField; 
            obj.variables2print{istre}.d_u          = reshape(vars.d_u',ndimField,[])';
            obj.variables2print{istre}.fext         = reshape(vars.fext',ndimField,[])';
            obj.variables2print{istre}.stress       = vars.stress;
            obj.variables2print{istre}.strain       = vars.strain;
            obj.variables2print{istre}.stress_fluct = vars.stress_fluct;
            obj.variables2print{istre}.strain_fluct = vars.strain_fluct;
        end
        
    end

    methods (Access = protected)

        function f = createVariablesToPrint(obj)
            f = obj.variables2print;
        end

        function t = createPrintType(obj)
           t = 'HomogenizedTensor';
        end

    end    

end