function surrogate = buildSidesSurrogate(problemNN,s)
model = problemNN;
if strcmp(s.output.type,'cholesky')
    model = CholeskySidesSurrogate(model,s.reference.trainFileName,s);
end
surrogate = SymmetricLatticeSurrogate(model,s.reference.trainFileName,s);
end