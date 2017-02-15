function out = getCoordVect(nPoints, samplDist)
    % A static method that returns an array of nPoints with a spacing of
    % sampleDist starting from -nPoints/2 to nPoints/2.
    out = -nPoints/2*samplDist:samplDist:(nPoints/2-1)*samplDist;
end