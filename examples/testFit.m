% Test fitting

close all
clear all

numRuns = 1;

for k = 1:numRuns
    angleInclination = pi/16;

    par.position = Length([0 0 0].*100,'nm');
    par.dipole = Dipole(angleInclination, 0);
    par.phaseMask = @(n) Vortex(n).cutInnerRing(0.1);
    par.stageDrift = LinearDrift(10,1,100,'nm');
    %par.defocus = Length(-500+1000*rand(), 'nm');

    psf = PSF(par);

    parEst.angleInclinationEstimate = angleInclination;
    parEst.angleAzimuthEstimate = 0;
    parEst.stageDrift = par.stageDrift;%BrownianMotion(10,30,'nm');
    fitResult = FitPSF(psf, parEst);
end

figure
plot(fitResult)