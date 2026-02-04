classdef CohesivePPRLaw < handle

    properties (Access = private)
        sigmaMax
        tauMax

        normalCharLength
        tangentialCharLength

        alpha
        beta       

        normalToughness
        tangentialToughness

        m
        n

        tauN
        tauT

        lambdaN
        lambdaT
    end

    methods (Access = public)

        function obj = CohesivePPRLaw(cParams)
            obj.init(cParams);
        end

        function t = evaluate(obj, disp)

            dn = disp(1);
            dt = disp(2);

            dnBar = dn / obj.normalCharLength;
            dtBar = dt / obj.tangentialCharLength;


            %% ???????????????????????????????????
            tn = 0; 
            tt = 0;  

            t = [tn; tt];
        end

        function D = derivative(obj, disp)

            dn = disp(1);
            dt = disp(2);

            dnBar = dn / obj.normalCharLength;
            dtBar = dt / obj.tangentialCharLength;


            % ??????????????????????????
            % dtn/ddn, dtn/ddt
            % dtt/ddn, dtt/ddt

            D = zeros(2,2); % placeholder
        end

    end


    methods (Access = private)

        function init(obj, cParams)
            obj.sigmaMax = cParams.sigmaMax;
            obj.tauMax   = cParams.tauMax;

            obj.normalCharLength     = cParams.normalCharLength;
            obj.tangentialCharLength = cParams.tangentialCharLength;

            obj.alpha = cParams.alpha;
            obj.beta  = cParams.beta;
    
            obj.normalToughness = cParams.normalToughness;
            obj.tangentialToughness = cParams.tangencialToughness;


            obj.m = obj.alpha*(obj.alpha-1)*obj.lambdaN^2 / (1-obj.alpha*obj.lambdaN^2);
            obj.n = obj.beta*(obj.beta-1)*obj.lambdaT^2 / (1-obj.beta*obj.lambdaT^2);


            obj.tauN = (-obj.normalToughness) ^ (obj.Macaulay(obj.normalToughness-obj.tangentialToughness)/ ...
                (obj.normalToughness-obj.tangentialToughness)) * (obj.beta/obj.m)^obj.m;
            obj.tauT = (-obj.tangentialToughness) ^ (obj.Macaulay(obj.tangentialToughness-obj.normalToughness)/ ...
                (obj.tangentialToughness-obj.normalToughness)) * (obj.beta/obj.n)^obj.n;

        end


        function m = Macaulay(a)
            if a<0 m=0;
            else m=a;
            end
        end

    end

end
