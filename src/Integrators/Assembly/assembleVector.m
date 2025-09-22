function F = assembleVector(Felem, f)
dofConnec = f.getDofConnec();
nDofs     = numel(f.fValues);
rowIdx    = dofConnec(:);
Felem = Felem';
F = sparse(rowIdx, 1, Felem(:), nDofs, 1);
end
