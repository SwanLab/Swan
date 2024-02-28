function dom = DDP(varargin)
    switch nargin
        case 2
            C = varargin{1};
            eps = varargin{2};
            s.operation = @(xV) evaluate2(C, eps, xV);
            dom = DomainFunction(s);
        case 3
            epsT = varargin{1};
            C    = varargin{2};
            eps  = varargin{3};
            s.operation = @(xV) evaluate3(epsT, C, eps, xV);
            dom = DomainFunction(s);
    end
end

function fVR = evaluate2(C, eps, xV)
    strn(:,1,:,:) = eps.evaluate(xV);
    Cv   = C.evaluate(xV);

    strs = pagemtimes(Cv,strn);
    % fVR = permute(strs, [1 3 4 2]);
    fVR = squeezeParticular(strs, 2);
end

function fVR = evaluate3(epsT, C, eps, xV)
    strnT(1,:,:,:) = epsT.evaluate(xV);
    Cv   = C.evaluate(xV);
    strn(:,1,:,:) = eps.evaluate(xV);

    strs = pagemtimes(strnT,Cv);
    E = pagemtimes(strs,strn);
    fVR = squeezeParticular(E,2);
end