function [psf, Ix, Iy] = getIntensitiesCamera(obj, fieldBFP)

    Ix = abs( obj.chirpZTransform.apply(obj, fieldBFP.x) ).^2; % intensity = |E_imagePlane|²
    Iy = abs( obj.chirpZTransform.apply(obj, fieldBFP.y) ).^2; % intensity = |E_imagePlane|²

    psf = Ix + Iy;

    % Account for pixel sensitivity
    psf = psf .* repmat(obj.pixelSensitivityMask, obj.nPixels);
    Ix = Ix .* repmat(obj.pixelSensitivityMask, obj.nPixels);
    Iy = Iy .* repmat(obj.pixelSensitivityMask, obj.nPixels);

    % Sum over block matrices (to return to desired pixelsize)
    [psf, Ix, Iy] = reduceOversampling(obj, psf, Ix, Iy);

    % Normalization
    totalIntensity = sum(sum(psf));
    psf = psf ./ totalIntensity * obj.nPhotons;
    Ix = Ix ./ totalIntensity * obj.nPhotons;
    Iy = Iy ./ totalIntensity * obj.nPhotons;

end