function par = readParameters(obj)

    par.nPixels = obj.nPixels;
    
    % Fluorophore
    par.dipole = obj.dipole;
    par.position = obj.position;
    par.nPhotons = obj.nPhotons;
    par.shotNoise = obj.shotNoise;
    par.reducedExcitation = obj.reducedExcitation;

    % Microscope setup
    par.wavelength = obj.wavelength;
    par.defocus = obj.defocus;
    par.astigmatism = obj.astigmatism;
    par.objectiveNA = obj.objectiveNA;
    par.objectiveFocalLength = obj.objectiveFocalLength;
    par.refractiveIndices = obj.refractiveIndices;
    par.heightIntermediateLayer = obj.heightIntermediateLayer;

    % Back focal plane
    par.phaseMask = obj.phaseMask;
    par.nDiscretizationBFP = obj.nDiscretizationBFP;

    % Camera
    par.pixelSize = obj.pixelSize;
    par.pixelSensitivityMask = obj.pixelSensitivityMask;
    par.backgroundNoise = obj.backgroundNoise;
end