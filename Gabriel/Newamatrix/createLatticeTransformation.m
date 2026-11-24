function transformation = createLatticeTransformation(settings)
type = 'Boost';
if isfield(settings,'cellTransformation')
    type = settings.cellTransformation;
end
switch type
    case 'Boost'
        transformation = BoostTransformation();
    case 'Shear'
        transformation = ShearTransformation();
    case 'SquareTurn'
        transformation = SquareTurnTransformation();
    case 'QuarterTurn'
        transformation = QuarterTurnTransformation();
    case 'Rotation'
        transformation = RotationTransformation();
    case 'TurnBoostStretch'
        transformation = TurnBoostStretchTransformation();
    otherwise
        error('createLatticeTransformation:UnknownType','Unknown cell transformation: %s',type);
end
end