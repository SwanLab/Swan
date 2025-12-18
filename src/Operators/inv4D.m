function Ainv = inv4D(A)
    Avoigt = convert2Voigt(A,'Material');
    AvoigtInv = inv(Avoigt);
    Ainv = covert2Tensor(AvoigtInv,'Material');
end