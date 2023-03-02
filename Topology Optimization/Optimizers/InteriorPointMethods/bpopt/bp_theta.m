function [th] = bp_theta(bp,x,s,bL,bU)
    th = sum(abs(bp_res(bp,x,s,bL,bU)));
