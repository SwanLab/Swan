function [L] = bp_lagrangian(bp,x,s,lam,bL,bU)

L = bp_obj(bp,x) + bp_res(bp,x,s,bL,bU)*lam';
