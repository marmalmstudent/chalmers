function [hFig, hSurf] = OctavePlot(X, xlab, Y, ylab, Z, zlab, n)
    hold('on');
    hFig = figure(n);
    hSurf = surf(X, Y, Z);
    xlabel(xlab);
    ylabel(ylab);
    zlabel(zlab);
    rotate3d('on'); # enable rotating
    view(30, 50); # azimuth, elevation
    shading('interp') # smooth surface, no grid
    colormap('jet')
end