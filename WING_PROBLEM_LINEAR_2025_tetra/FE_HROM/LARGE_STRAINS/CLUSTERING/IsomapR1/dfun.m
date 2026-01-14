function out = dfun(in)

global X

out = L2_distance(X,X(:,in)); 
out = out'; 
