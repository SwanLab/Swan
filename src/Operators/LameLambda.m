function dom = LameLambda(E,nu,N)
    kappa = E./(N*(1-(N-1)*nu));
    mu    = E./(2*(1+nu));
    dom   = kappa - (2/N)*mu;
end