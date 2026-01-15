classdef TrainedRVE < handle

    properties (GetAccess = public , SetAccess = private)
        ndimf
        Kcoarse
        Udef
        Urb
        U
        DOFl
        DOFr
        PhiDef
        PhiRb
        Grb
        dKcoarse
        dUrb
        dUdef
        dU
    end

    properties (Access = private)


    end

    properties (Access = private)

    end

    methods (Access = public)

        function obj = TrainedRVE(filename)
            obj.init(filename)

        end

    end

    methods (Access = private)
        
        function init(obj,filename)
            load(filename);
            obj.Kcoarse = EIFEoper.Kcoarse;
            if isfield(EIFEoper,'Udef')
                obj.Udef    = EIFEoper.Udef;
                obj.Urb     = EIFEoper.Urb;
            else
                obj.Udef    = EIFEoper.RECONSTRUCTION.DEF_DISP.BASIS*...
                              EIFEoper.RECONSTRUCTION.DEF_DISP.coeff;
                obj.Urb     = EIFEoper.RECONSTRUCTION.RB_DISP.BASIS*...
                              EIFEoper.RECONSTRUCTION.RB_DISP.coeff;
            end
            
            obj.ndimf   = 2;
            %             obj.PhiDef  = EIFEoper.RECONSTRUCTION.DEF_DISP.BASIS;
            %             obj.PhiRb   = EIFEoper.RECONSTRUCTION.RB_DISP.BASIS;
            %             obj.Grb     = EIFEoper.Grb;
        end

    end

end