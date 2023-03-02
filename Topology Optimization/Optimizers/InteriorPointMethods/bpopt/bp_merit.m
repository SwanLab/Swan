function [me] = bp_merit(bp,x,xL,xU,s,bL,bU)
    ph = bp_phi(bp,x,xL,xU,s,bL,bU);
    r = bp_res(bp,x,s,bL,bU);
    me = ph + bp.nu*sum(abs(r));