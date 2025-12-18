function Ainv = inv4D(A)
    Avoigt = convert2Voigt(A,'Material');
    AvoigtInv = inv(Avoigt);
    Ainv = convert2Tensor(AvoigtInv,'Material');
end