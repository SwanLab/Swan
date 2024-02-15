function mesh = UnitHexaMesh(nx, ny, nz)
    length = 1;
    height = 1;
    width  = 1;
    mesh = HexaMesh(length, height, width, nx, ny, nz);
end