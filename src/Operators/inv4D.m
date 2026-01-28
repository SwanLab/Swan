function Ainv = inv4D(A)
    Avoigt = convert2Voigt(A,'Constitutive');
    AvoigtInv = inv(Avoigt);
    Ainv = convert2Tensor(AvoigtInv,'Compliance');
end