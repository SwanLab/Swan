classdef PostprocessDensityInGaussPoints < Postprocess_TopOpt
    
    
    properties (Access = private)
        FieldName     = 'Density';
        ComponentName = 'Density'
    end
    
    
    methods (Access = public)
        
        function obj = PostprocessDensityInGaussPoints(quadrature)
            obj.SetQuadratureInfo(quadrature)
        end
        
    end

    methods (Access = protected)    
        
        function PrintResultsOld(obj)
            obj.ngaus = size(obj.Field,2);
            obj.PrintGaussPointsHeader()
            obj.PrintScalar(obj.FieldName,obj.ComponentName,'Elastic Problem','Scalar','OnGaussPoints',obj.gauss_points_name,obj.Field,obj.Iter);
        end
        
        function printResults(obj)
            iD = obj.fid_res;
            fN = obj.file_name;
            nS = obj.nsteps;
            gD = obj.gauss_points_name;
            eT = obj.etype;
            pT = obj.ptype;
            nG = obj.ngaus;
            nD = obj.ndim;
            pG = obj.posgp;
            rS = obj.Field;
            iT = obj.Iter;
            DensityGaussResultsPrinter(iD,fN,nS,gD,eT,pT,nG,nD,pG,rS,iT);
        end
       
    end    
    
    methods (Access = private)
        function SetQuadratureInfo(obj,quadrature)
           obj.ngaus = quadrature.ngaus;
           obj.posgp = quadrature.posgp';        
        end
    end
    
end

