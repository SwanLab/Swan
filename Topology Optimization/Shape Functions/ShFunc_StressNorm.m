classdef ShFunc_StressNorm < ShFunWithElasticPdes
    
    properties (Access = private)
        adjointProb
        fAdjoint
    end
    
    methods (Access = public)
        
        function obj = ShFunc_StressNorm(cParams)
            cParams.filterParams.quadratureOrder = 'LINEAR';
            obj.init(cParams);
            fileName = cParams.femSettings.fileName;            
            obj.createEquilibriumProblem(fileName);
            obj.createAdjointProblem(fileName);
            obj.createOrientationUpdater();                        
        end
        
        function f = getPhysicalProblems(obj)
            f{1} = obj.physicalProblem;
            f{2} = obj.adjointProb;
        end
        
    end
    
    methods (Access = protected)
        
        function computeFunctionValue(obj)
            phy = obj.physicalProblem;
            dvolum  = phy.geometry.dvolu;
            stress = phy.variables.stress;
            P = obj.homogenizedVariablesComputer.P;
            ngaus = phy.element.quadrature.ngaus;
            nstre = size(stress,2);
            c = 0;
            for igaus = 1:ngaus
                stressG = squeeze(stress(igaus,:,:));
                dV = dvolum(:,igaus);
                for istre = 1:nstre
                    sI(:,1) = stressG(istre,:);
                    for jstre = 1:nstre
                        sJ(:,1) = stressG(jstre,:);
                        Pij = squeeze(P(istre,jstre,:));
                        csum = sI.*Pij.*sJ.*dV;
                        c = c + sum(csum);
                    end
                end
            end
            obj.value = c;
        end
        
        function updateHomogenizedMaterialProperties(obj)
            obj.filterDesignVariable();
            obj.homogenizedVariablesComputer.computeCtensor(obj.regDesignVariable);
            obj.homogenizedVariablesComputer.computePtensor(obj.regDesignVariable);
        end
        
        function computeGradient(obj)
            g1 = obj.computeFirstTerm();
            g2 = obj.computeSecondTerm();
            g3 = obj.computeThirdTerm();
            g = g1 + g2 + g3;            
            gf = zeros(size(obj.Msmooth,1),obj.nVariables);
            for ivar = 1:obj.nVariables
                gs = squeeze(g(:,:,ivar));
                gf(:,ivar) = obj.filter.getP1fromP0(gs);
            end
            g = obj.Msmooth*gf;
            obj.gradient = g(:);
        end
        
        function g = computeFirstTerm(obj)
            nelem = obj.physicalProblem.geometry.interpolation.nelem;
            ngaus = obj.physicalProblem.element.quadrature.ngaus;
            nstre = obj.physicalProblem.element.getNstre();
            eu = obj.physicalProblem.variables.strain;
            dC = obj.homogenizedVariablesComputer.dC;
            C  = obj.homogenizedVariablesComputer.C;
            P  = obj.homogenizedVariablesComputer.P;
            g = zeros(nelem,ngaus,obj.nVariables);
            for igaus = 1:ngaus
                for istre = 1:nstre
                    ei   = squeeze(eu(igaus,istre,:));
                    for jstre = 1:nstre
                        for ivar = 1:obj.nVariables
                            dCij_iv = squeeze(dC(istre,jstre,ivar,:));
                            for kstre = 1:nstre
                                Pjk  = squeeze(P(jstre,kstre,:));
                                for lstre = 1:nstre
                                    Ckl  = squeeze(C(kstre,lstre,:));
                                    el   = squeeze(eu(igaus,lstre,:));
                                    g_iv  = squeeze(g(:,igaus,ivar));
                                    g(:,igaus,ivar) = g_iv + 2*ei.*dCij_iv.*Pjk.*Ckl.*el;
                                end
                            end
                        end
                    end
                end
            end
        end
        
        function g = computeSecondTerm(obj)
            nelem = obj.physicalProblem.geometry.interpolation.nelem;
            ngaus = obj.physicalProblem.element.quadrature.ngaus;
            nstre = obj.physicalProblem.element.getNstre();
            eu = obj.physicalProblem.variables.strain;
            C  = obj.homogenizedVariablesComputer.C;
            dP  = obj.homogenizedVariablesComputer.dP;
            g = zeros(nelem,ngaus,obj.nVariables);
            for igaus = 1:ngaus
                for istre = 1:nstre
                    ei   = squeeze(eu(igaus,istre,:));
                    for jstre = 1:nstre
                        for ivar = 1:obj.nVariables
                            Cij = squeeze(C(istre,jstre,:));
                            for kstre = 1:nstre
                                dPjk_iv  = squeeze(dP(jstre,kstre,ivar,:));
                                for lstre = 1:nstre
                                    Ckl  = squeeze(C(kstre,lstre,:));
                                    el   = squeeze(eu(igaus,lstre,:));
                                    g_iv  = squeeze(g(:,igaus,ivar));
                                    g(:,igaus,ivar) = g_iv + ei.*Cij.*dPjk_iv.*Ckl.*el;
                                end
                            end
                        end
                    end
                end
            end
        end
        
        function g = computeThirdTerm(obj)
            nelem = obj.physicalProblem.geometry.interpolation.nelem;
            ngaus = obj.physicalProblem.element.quadrature.ngaus;
            nstre = obj.physicalProblem.element.getNstre();
            eu = obj.physicalProblem.variables.strain;
            dC = obj.homogenizedVariablesComputer.dC;
            ep = obj.adjointProb.variables.strain;
            g = zeros(nelem,ngaus,obj.nVariables);
            for igaus = 1:ngaus
                for istre = 1:nstre
                    ei   = squeeze(eu(igaus,istre,:));
                    for jstre = 1:nstre
                        ej   = squeeze(ep(igaus,jstre,:));
                        for ivar = 1:obj.nVariables
                            dCij_iv = squeeze(dC(istre,jstre,ivar,:));                            
                            g_iv  = squeeze(g(:,igaus,ivar));
                            g(:,igaus,ivar) = g_iv + ei.*dCij_iv.*ej;
                        end
                    end                    
                end
            end
        end
        
        function solvePDEs(obj)
            obj.physicalProblem.setC(obj.homogenizedVariablesComputer.C);
            obj.physicalProblem.computeVariables();
            obj.computeFadjoint();
            obj.adjointProb.setC(obj.homogenizedVariablesComputer.C);
            obj.adjointProb.computeVariablesWithBodyForces(obj.fAdjoint);
        end
        
        function updateGradient(obj)
            
            
        end
    end
    
    methods (Access = private)
        
        function createAdjointProblem(obj,fileName)
            obj.adjointProb = FEM.create(fileName);
            obj.adjointProb.preProcess;
        end
        
        function computeFadjoint(obj)
            phy    = obj.physicalProblem;
            dvolum = phy.geometry.dvolu;
            strain = phy.variables.strain;
            P = obj.homogenizedVariablesComputer.P;
            C = obj.homogenizedVariablesComputer.C;
            ngaus = phy.element.quadrature.ngaus;
            nstre = size(strain,2);
            nelem = size(strain,3);
            nnode = phy.element.nnode;
            nunkn = phy.element.dof.nunkn;
            
            CPC = zeros(size(C));
            for istre = 1:nstre
                for jstre = 1:nstre
                    for kstre = 1:nstre
                        for lstre = 1:nstre
                            CPC(istre,jstre,:) = CPC(istre,jstre,:) + C(istre,kstre,:).*P(kstre,lstre,:).*C(lstre,jstre,:);
                        end
                    end
                end
            end
            
            eforce = zeros(nunkn*nnode,ngaus,nelem);
            for igaus = 1:ngaus
                Bmat = phy.element.computeB(igaus);
                strainG = squeeze(strain(igaus,:,:));
                dV(:,1) = dvolum(:,igaus);
                for istre = 1:nstre
                    ei(:,1) = squeeze(strainG(istre,:));
                    for jstre = 1:nstre
                        CPCij = squeeze(CPC(istre,jstre,:));
                        for iv = 1:nnode*nunkn
                            Bj_iv = squeeze(Bmat(jstre,iv,:));
                            int = -2*ei.*CPCij.*Bj_iv;
                            Fiv = squeeze(eforce(iv,igaus,:));
                            eforce(iv,igaus,:) = Fiv + int.*dV;
                        end
                    end
                end
            end
            Fvol = phy.element.AssembleVector({eforce});
            obj.fAdjoint = Fvol;
        end
        
    end
    
end
