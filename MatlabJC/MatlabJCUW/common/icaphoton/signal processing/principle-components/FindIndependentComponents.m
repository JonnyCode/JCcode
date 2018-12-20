function ReturnEigVec = FindIndependentComponents(EigVec, slope)

theta = atan(slope)

ReturnEigVec = EigVec;
ReturnEigVec(:, 1) = EigVec(:, 1) .* cos(theta) + EigVec(:, 2) .* sin(theta);
ReturnEigVec(:, 2) = -EigVec(:, 1) .* sin(theta) + EigVec(:, 2) .* cos(theta);
