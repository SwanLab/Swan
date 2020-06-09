classdef VademecumReader < handle
    
    properties (Access = public)
        qV
        phiV
        xiV
        rhoV
        mxV
        myV
    end
    
    properties (Access = private)
        iter
        nMx
        nMy
        nPhi
    end
    
    properties (Access = private)
        compressedFileName
        vademecum
    end
    
    methods (Access = public)
        
        function obj = VademecumReader()
            obj.init();
            obj.loadCompressedVademecum();
            obj.transformDataInMatrixToVectorForm();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.compressedFileName = 'OptimalSuperEllipses';
        end
        
        function loadCompressedVademecum(obj)
            d = load(obj.compressedFileName);
            obj.vademecum = d;
        end
        
        function transformDataInMatrixToVectorForm(obj)
            obj.computeNmxNmyNphi();
            obj.initVademecumVariables();
            for iMx = 1:obj.nMx
                for iMy = 1:obj.nMy
                    obj.computeVademecumVariables(iMx,iMy);
                end
            end
        end
        
        function computeNmxNmyNphi(obj)
            v = obj.vademecum;
            obj.nMx  = size(v.q,1);
            obj.nMy  = size(v.q,2);
            obj.nPhi = size(v.q,3);
        end
        
        function initVademecumVariables(obj)
            nP = obj.nPhi;
            nT   = obj.nMx*obj.nMy;
            obj.qV   = zeros(nT,nP);
            obj.phiV = zeros(nT,nP);
            obj.mxV  = zeros(nT,nP);
            obj.myV  = zeros(nT,nP);
            obj.xiV  = zeros(nT,1);
            obj.rhoV = zeros(nT,1);            
        end
        
        function computeVademecumVariables(obj,iMx,iMy)
            obj.iter = (iMx-1)*obj.nMy + iMy;
            obj.computeQ(iMx,iMy);
            obj.computePhi(iMx,iMy);
            obj.computeXi(iMx,iMy);
            obj.computeRho(iMx,iMy);
            obj.computeMx();
            obj.computeMy();
        end
        
        function computeQ(obj,iMx,iMy)
            v = obj.vademecum;
            q = v.q(iMx,iMy,:);
            obj.qV(obj.iter,:) = q;
        end
        
        function computePhi(obj,iMx,iMy)
            v = obj.vademecum;
            phi = v.phi(iMx,iMy,:);
            obj.phiV(obj.iter,:) = phi;
        end
        
        function computeXi(obj,iMx,iMy)
            v = obj.vademecum;
            xi = v.xi(iMx,iMy,:);
            obj.xiV(obj.iter,1) = unique(xi);
        end
        
        function computeRho(obj,iMx,iMy)
            v = obj.vademecum;
            rho =  v.rho(iMx,iMy,:);
            obj.rhoV(obj.iter,1) = unique(rho);
        end
        
        function computeMx(obj)
            xi  = obj.xiV(obj.iter,1);
            rho = obj.rhoV(obj.iter,1);
            q   = obj.qV(obj.iter,:);
            mx  = SuperEllipseParamsRelator.mx(xi,rho,q);
            obj.mxV(obj.iter,:) = mx;
        end
        
        function computeMy(obj)
            xi  = obj.xiV(obj.iter,1);
            rho = obj.rhoV(obj.iter,1);
            q   = obj.qV(obj.iter,:);
            my  = SuperEllipseParamsRelator.my(xi,rho,q);
            obj.myV(obj.iter,:) = my;
        end
        
    end
    
end