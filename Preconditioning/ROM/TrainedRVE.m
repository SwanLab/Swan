classdef TrainedRVE < handle

    properties (GetAccess = public , SetAccess = private)
        ndimf
        Kcoarse
        Udef
        Urb
        U
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

        function init(obj,data)
                load(data);
                obj.Kcoarse = EIFEoper.Kcoarse;
                if isfield(EIFEoper,'Udef')
                    obj.Udef    = EIFEoper.Udef;
                    obj.Urb     = EIFEoper.Urb;
                elseif isfield(EIFEoper, 'RECONSTRUCTION')
                    obj.Udef    = EIFEoper.RECONSTRUCTION.DEF_DISP.BASIS*...
                                  EIFEoper.RECONSTRUCTION.DEF_DISP.coeff;
                    obj.Urb     = EIFEoper.RECONSTRUCTION.RB_DISP.BASIS*...
                                  EIFEoper.RECONSTRUCTION.RB_DISP.coeff;
                end
                
                if isfield(EIFEoper,'U')
                    obj.U = EIFEoper.U;
                end
                if isfield(EIFEoper,'dKcoarse')
                    obj.dKcoarse = EIFEoper.dKcoarse;
                end
                if isfield(EIFEoper,'dUdef')
                    obj.dUdef = EIFEoper.dUdef;
                end
                if isfield(EIFEoper,'dUrb')
                    obj.dUrb = EIFEoper.dUrb;
                end
                if isfield(EIFEoper,'dU')
                    obj.dU = EIFEoper.dU;
                end

            obj.ndimf   = 2;
            %             obj.PhiDef  = EIFEoper.RECONSTRUCTION.DEF_DISP.BASIS;
            %             obj.PhiRb   = EIFEoper.RECONSTRUCTION.RB_DISP.BASIS;
            %             obj.Grb     = EIFEoper.Grb;
        end

    end

end