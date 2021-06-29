classdef AmplificatorComponentsCalculator < handle
    
    properties (SetAccess = private, GetAccess = public)
        Phomog
        monom        
    end
    
    properties (Access = private)
        integrationDB
        Ptensor
        nQuadStre
        monomCoef
        nMonom
        pExp     
        PcoefHomog    
        tstress
        Chomog
    end
    
    methods (Access = public)
        
        function obj = AmplificatorComponentsCalculator(d)
            obj.init(d)
        end
        
        function compute(obj)
            obj.computeQuadStressComponents();
            obj.computeMonomials();
            obj.computeMonomialsCoeficients()
            obj.integratePtensor();
            obj.computePtensorWithCoef();
        end
        
    end
    
    methods (Access = private)
                                                                                                                                                                                                                    
        function init(obj,d)
            obj.integrationDB.nstre = d.nstre;
            obj.integrationDB.V = d.V;
            obj.integrationDB.dV = d.dV;
            obj.integrationDB.ngaus = d.ngaus;
            obj.Chomog  = d.Ch;
            obj.tstress = d.tstress;
            obj.pExp  = d.pNorm/2;
        end
        
        function computeMonomials(obj)
            p = obj.pExp;
            n = obj.nQuadStre;
            alpha = multinomial_expand(p,n);
            obj.monom  = alpha;
            obj.nMonom = size(alpha,1);
        end
        
        function computeMonomialsCoeficients(obj)
            a = obj.monom;
            p = obj.pExp;            
            ncoefs = size(a,1);
            coef = zeros(ncoefs,1);
            for icoef = 1:ncoefs
                coef(icoef,1) = multcoef(p,a(icoef,:));
            end
            obj.monomCoef = coef;
        end
        
        function computeQuadStressComponents(obj)
            nstre = obj.integrationDB.nstre;
            nQ = nstre*(nstre+1)/2;
            obj.nQuadStre = nQ;
        end
        
        function integratePtensor(obj)
            P2v = obj.computeP2inVector();
            V     = obj.integrationDB.V;
            dV    = obj.integrationDB.dV;
            ngaus = obj.integrationDB.ngaus;
            nelem = size(P2v,3);
            prodTerm  = zeros(ngaus,nelem);            
            P = zeros(obj.nMonom,1);
            for t = 1:obj.nMonom
                integrand = ones(ngaus,nelem);
                for k = 1:obj.nQuadStre                  
                   prodTerm(1:ngaus,:) = squeeze(P2v(k,:,:));
                   alpha = obj.monom(t,k);
                   integrand = integrand.*(prodTerm.^alpha);
                end
                int = integrand.*dV;
                P(t) = 1/V*sum(int(:));               
            end
            obj.Phomog = P;
        end
        
        function P2 = computeP2(obj)            
            Pt    = obj.computePtensor();
            ngaus = size(Pt,2);
            nstre = size(Pt,1);
            nelem = size(Pt,4);
            P2 = zeros(nstre,nstre,ngaus,nelem);            
            for istre = 1:nstre
                for jstre = 1:nstre
                    P2ij = zeros(ngaus,nelem);
                    for kstre = 1:nstre
                         Pki(1:ngaus,:) = squeeze(Pt(kstre,:,istre,:));
                         Pkj(1:ngaus,:) = squeeze(Pt(kstre,:,jstre,:));
                         if kstre == 3
                             f = 2;
                         else
                             f = 1;
                         end
                         P2ij = P2ij + f*Pki.*Pkj;
                    end
                    P2(istre,jstre,1:ngaus,:) = P2ij;
                end
            end                        
        end
        
        function P2v = computeP2inVector(obj)
            P2 = obj.computeP2();
            nstre = size(P2,1);
            ngaus = size(P2,3);
            nelem = size(P2,4);
            nQuads = nstre*(nstre+1)/2;
            P2v = zeros(nQuads,ngaus,nelem);
            for s = 1:nQuads
                   [istre,jstre]    = obj.indexVoigt2tensor(s);
                   factor           = obj.computeVoigtFactor(s,nstre);
                   P2v(s,1:ngaus,:) = factor*squeeze(P2(istre,jstre,:,:));       
             end                         
        end
        
        function computePtensorWithCoef(obj)
            obj.PcoefHomog = obj.Phomog.*obj.monomCoef;
        end
        
        function Pt = computePtensor(obj)
            tstres = obj.tstress;
            Shomog = obj.computeShomog();
            nstre  = obj.integrationDB.nstre;
            Pt = zeros(size(tstres));
            for istre = 1:nstre
                for jstre = 1:nstre
                    for kstre = 1:nstre
                        Pt(istre,:,jstre,:) = Pt(istre,:,jstre,:) + tstres(istre,:,kstre,:)*Shomog(kstre,jstre);
                    end
                end
            end
        end           
        
        function Shomog = computeShomog(obj)
            Ct = StiffnessPlaneStressVoigtTensor;
            Ct.setValue(obj.Chomog);
            St = Inverter.invert(Ct);
            Shomog = St.getValue();            
        end        
        
    end
    
    methods (Access = private, Static)
        
        function f = computeVoigtFactor(k,nstre)
            if nstre == 3
                if k <= 3
                    f = 1;
                else
                    f = 2;
                end
            elseif nstre == 6
                if k <= 3
                    f = 1;
                else
                    f = 2;
                end
            end
        end
        
        function [i,j] = indexVoigt2tensor(k)
            T = [1 1;
                2 2;
                3 3;
                2 3;
                1 3;
                1 2];
            i = T(k,1);
            j = T(k,2);
        end
       
    end
    
    
end



