function [E2, L1, L2, uvect, vvect, umat, vmat]=TSM(E1, a, b, L, lbda_0, n)
    % This function implements the conventional Huygens-Fresnel method (TSM) for
    % free space propagation.
    N = size(E1, 1);  # Matrix size (one side of the square matrix E1)
    lambda_medium = lbda_0/n;
    k = 2*pi/lambda_medium;

    # Plane 1 coordinates:
    xvect = getCoordVect(N, a);
    yvect = getCoordVect(N, a);
    [xmat, ymat] = meshgrid(xvect, yvect);

    # Plane 2 coordinates
    uvect = getCoordVect(N, b);
    vvect = getCoordVect(N, b);
    [umat, vmat] = meshgrid(uvect, vvect);

    # Distance to dummy plane
    L1 = L*a/(a-b);  # plane 1 to dummy plane
    L2 = L*b/(a-b);  # plane 2 to dummy plane
    # Sampling distance in Dummy plane
    c = lambda_medium/(N*a)*L1;

    # Dummy plane coordinates
    xiVect = getCoordVect(N, c);
    etaVect = getCoordVect(N, c);
    [xiMat, etaMat] = meshgrid(xiVect, etaVect);

    # Calculating the field in Plane 2
    fctr0 = a^2 / b^2 * L2/L1 * exp(1j*k*(L1-L2));
    fctr1 = exp(-1j*k *
                (umat.^2 + vmat.^2) / (2*L2));
    prefactor = exp(
        1j*k * (xiMat.^2 + etaMat.^2) / 2 * (1/L1 - 1/L2));
    fft_part = fft2c(
        E1.*exp(1j*k * (xmat.^2 + ymat.^2) / (2*L1)));
    fctr2 = ifft2c(prefactor .* fft_part);
    E2 = fctr0.*fctr1.*fctr2;
    % return E2, L1, L2, uvect, vvect, umat, vmat
end