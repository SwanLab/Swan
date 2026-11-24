function surrogate = buildLatticeSurrogate(problemNN,s)
model = problemNN;
if strcmp(s.output.type,'cholesky')
    model = CholeskyLatticeSurrogate(model,s.reference.trainFileName,s);
end
surrogate = SymmetricLatticeSurrogate(model,s.reference.trainFileName,s);
end