classdef LHSintergrator_StiffnessElasticStoredB < LHSintegrator

    properties (Access = private)
        geometry
        Btot
    end

    methods (Access = public)
        
        function obj = LHSintergrator_StiffnessElasticStoredB(cParams)
            obj.init(cParams)
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
            obj.computeB();
        end

        function LHS = compute(obj)
            CmatTot = obj.assemblyCmat();
            LHS = obj.computeStiffness(CmatTot);
        end
        
        function setMaterialC(obj, Cmat)
            obj.material.C = Cmat;
        end
    end
    
   methods (Access = protected)
        
        function lhs = computeElementalLHS(obj)
        end

        function lhs = assembleMatrix(obj)
        end
        
   end
    
   methods (Access = private)

       function createGeometry(obj)
           s.mesh = obj.mesh;
           g      = Geometry.create(s);
           quad   = obj.quadrature;
           interp = obj.interpolation;
           g.computeGeometry(quad,interp); 
           obj.geometry = g;
       end

       function computeB(obj)
           s.dim          = obj.dim;
           s.geometry     = obj.geometry;
           s.globalConnec = obj.globalConnec;
           BMC  = BMatrixComputer(s);
           obj.Btot = BMC.compute();
       end

       function CmatTot = assemblyCmat(obj)
           nstre = obj.dim.nstre;
           ngaus = obj.dim.ngaus;
           ntot  = obj.dim.nt;
           Cmat    = obj.material.C;
           CmatTot = sparse(ntot,ntot);
           dvol = obj.geometry.dvolu;
           for istre = 1:nstre
               for jstre = 1:nstre
                   for igaus = 1:ngaus

                       posI = (istre)+(nstre)*(igaus-1) : ngaus*nstre : ntot;
                       posJ = (jstre)+(nstre)*(igaus-1) : ngaus*nstre : ntot;
                       
                       Ct = squeeze(Cmat(istre,jstre,:,igaus)).*dvol(:,igaus);
                       CmatTot = CmatTot + sparse(posI,posJ,Ct,ntot,ntot);
                   end
               end
           end
       end

       function K = computeStiffness(obj,CmatTot)
           B  = obj.Btot;
           CB = CmatTot*B;
           K  = B'*CB;
       end

   end
    
end