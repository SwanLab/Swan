function MultinomialTest

p = 2;
n = 6; %s1s1,s2s2,s12s12,s12s2,s1s12,s1s2

alpha = multinomial_expand(p,n); 
ncoefs = size(alpha,1);

coef = zeros(ncoefs,1);
for icoef = 1:ncoefs
coef(icoef,1) = multcoef(p,alpha(icoef,:));
end

[alpha,coef]

end