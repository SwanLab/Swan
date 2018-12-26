classdef Postprocess_TopOpt_density < Postprocess_TopOpt
    
    properties (Access = private)
        FieldName = 'Density';
    end
            
   methods (Access = protected)
        
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
            DensityResultsPrinter(iD,fN,nS,gD,eT,pT,nG,nD,pG,rS,iT);
        end
   end
    
end
