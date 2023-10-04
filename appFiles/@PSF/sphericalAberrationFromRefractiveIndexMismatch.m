function [Defocus, SA, a] = sphericalAberrationFromRefractiveIndexMismatch(obj, removeDefocus)

% Outputs spherical aberration arising from a RI-mismatch, for z=1m inside
% the material of RI2. Output can be scaled to the existing z-value
% returned defocus phase is defined for medium n2.
% According to paper: Jesacher et al. 2010. Opex
% The "defocus-term" is removed analytically.
%
% Inputs:
% n1 ... refractive index of immersion medium
% n2 ... refractive index of target material
% NA ... numerical aperture of objective
% lambda_0 ... vacuum wavelength
% a  ... coefficient of "Defocus" contained in SA
% rd ... 0 or 1, if set to 1, defocus is removed

N = obj.nDiscretizationBFP;
NA = obj.objectiveNA;
n2 = obj.refractiveIndices(3);
k = 2*pi / obj.wavelength.inMeter;

par.nGrid = N;
par.spacing = 2/(N-1);
par.mode = 'FFT';
mask = Mask(par);
[normR, ~] = getPolarCoordinates(mask);

% Mean value of spherical defocus
MW = @(RI) 2/3*k*(-(-1 + RI^2/NA^2)^(3/2)+(RI^2/NA^2)^(3/2))*NA;
% Defocus function
Def = @(RI) real(k*NA*sqrt(RI^2/NA^2-normR.^2)-MW(RI)).*mask.values; %spherical defocus in medium with refractive index RI (for 1m of defocus)
Defocus = Def(n2);

if nargout > 1
    n1 = obj.refractiveIndices(3);
    SA = -(Def(n2)-Def(n1)); %spherical aberration phase (complete, i.e including defocus) (for 1m of defocus), (see paper: Jesacher & Booth, Opex 2010)

    a_spher = -1/72*k^2*pi*(72*n2^2-36*NA^2+(32*(-(-1 + n2^2/NA^2)^(3/2) + n2^3/NA^3)*(n1^3-n1^2*sqrt((n1 - NA)*(n1 + NA))+NA^2*sqrt((n1 - NA)*(n1 + NA))))/NA - (32*(n2^3 - n2^2*sqrt((n2 - NA)*(n2 + NA))+NA^2*sqrt((n2 - NA)*(n2 + NA)))^2)/NA^4+(9/NA^2)*(2*(-n1^3*n2 - n1*n2^3 + n1^2*sqrt((n1 - NA)*(n2 - NA)*(n1 + NA)*(n2 + NA))+sqrt((n1 - NA)*(n2 - NA)*(n1 + NA)*(n2 + NA))*(n2^2-2*NA^2)) - (n1^2 - n2^2)^2*(log((n1 - n2)^2)-log(n1^2 + n2^2 - 2*(NA^2 + sqrt((n1 - NA)*(n2 - NA)*(n1 + NA)*(n2 + NA)))))));
    def_norm = -((k^2*(16*n2^6 - 24*n2^4*NA^2 + 6*n2^2*NA^4 + NA^6 - 16*n2^5*sqrt(n2^2 - NA^2) + 16*n2^3*NA^2*sqrt(n2^2 - NA^2))*pi)/(18*NA^4));
    a = a_spher/def_norm; %coefficient of Defocus contained in SA (without normalization)

    if removeDefocus
        SA = SA - a_spher*Def(n2)/def_norm;
    end
end

end