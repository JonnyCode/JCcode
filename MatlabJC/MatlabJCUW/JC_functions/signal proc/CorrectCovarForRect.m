function [CovarRectDist, MeansRectDist] = CorrectCovarForRect(OrigCovar, OrigMeans, RectThresholds)


% project to eigenvector space
[EigVec, EigVal] = eig(OrigCovar);
for dim = 1:length(OrigMeans)
    if (-min(EigVec(dim, :)) > max(EigVec(dim, :)))
        EigVec(dim, :) = -EigVec(dim, :);
    end
end
ProjectedMeans = OrigMeans * EigVec;
ProjectedThresholds = RectThresholds * EigVec;
SampSize = 0.01;

% compute rectified moments
NewCovar = zeros(length(OrigMeans), length(OrigMeans));
for dim = 1:length(OrigMeans)
    [NewProjMeans(dim), NewProjSD(dim)] = RectGaussMoments(ProjectedMeans(dim), sqrt(EigVal(dim, dim)), ProjectedThresholds, SampSize);
    NewCovar(dim, dim) = NewProjSD(dim)^2;
end

% project back to original space
CovarRectDist = inv(EigVec) * NewCovar * EigVec;
MeansRectDist = NewProjMeans * inv(EigVec);
