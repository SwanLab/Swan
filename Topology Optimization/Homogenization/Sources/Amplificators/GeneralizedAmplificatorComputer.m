classdef GeneralizedAmplificatorComputer < handle
    
    properties (SetAccess = private, GetAccess = public)
        Phomog
    end
    
    properties (Access = private)
        integrationDB
        Ptensor
        pNorm
        nQuadStre
        monom
        monomCoef
        nMonom
        pExp     
        PcoefHomog        
    end
    
    methods (Access = public)
        
        function obj = GeneralizedAmplificatorComputer(d)
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
            obj.integrationDB.ngaus  = d.ngaus;
            obj.Ptensor = d.Ptensor;
            obj.pNorm = 4;
            obj.pExp  = obj.pNorm/2;
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
            n = 0;
            for i = 1:obj.integrationDB.nstre
                n = n + i;
            end
            obj.nQuadStre = n;
        end
        
        function integratePtensor(obj)
            nstre = obj.integrationDB.nstre;
            V     = obj.integrationDB.V;
            dV    = obj.integrationDB.dV;
            ngaus = obj.integrationDB.ngaus;
            Pt    = obj.Ptensor;
            
            P = zeros(obj.nMonom,1);
            for t = 1:obj.nMonom
                integrand = 0;                
                for igaus = 1:ngaus
                    integrandG = ones(size(Pt,4),1);
                    for s = 1:obj.nQuadStre
                        [istre,jstre] = obj.indexVoigt2tensor(s);
                        prodTerm = zeros(size(Pt,4),1);
                        for kstre = 1:nstre
                            Pki = squeeze(Pt(kstre,igaus,istre,:));
                            Pkj = squeeze(Pt(kstre,igaus,jstre,:));
                            factor = obj.computeVoigtFactor(kstre,nstre);
                            prodTerm = prodTerm + factor*(Pki.*Pkj);                                                        
                        end
                        alpha = obj.monom(t,s);
                        integrandG = integrandG.*(prodTerm.^alpha);
                    end
                   integrand = integrand + 1/V*integrandG'*dV(:,igaus);
                end
                P(t) = integrand;                
            end
            obj.Phomog = P;
        end
        
        function computePtensorWithCoef(obj)
            obj.PcoefHomog = obj.Phomog.*obj.monomCoef;
        end
        
    end
    
    methods (Access = private, Static)
        
        function f = computeVoigtFactor(k,nstre)
            if nstre == 3
                if k < 3
                    f = 1;
                else
                    f = 2;
                end
            elseif nstre ==6
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



