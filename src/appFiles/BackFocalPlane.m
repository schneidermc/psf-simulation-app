classdef BackFocalPlane

    % Calculates the BFP-fields (Ex,Ey) of an arbitrary oriented dipole in
    % the focal point of an objective based on the following paper:
    % Axelrod, J. of Microscopy 2012
    %
    % Details:
    % - An intermediate layer between coverglass and specimen can be assumed
    % - Includes also near-field effects (SAF emission)

    properties
        nGrid (1,1) {mustBeInteger, mustBePositive} = 129
        unitKSpace (1,1) {mustBePositive} = 1
        electricField
    end

    methods
        function obj = BackFocalPlane(psfObj)
            obj.nGrid = psfObj.nDiscretizationBFP;
            obj.unitKSpace = psfObj.unitKSpace; % k-space unit in BFP
            obj.electricField = calculateElectricField(obj,psfObj);
        end

        function E_BFP = calculateElectricField(obj,psfObj)

            RI = psfObj.refractiveIndices;
            dipole = psfObj.dipole;
            hIntermediate = psfObj.heightIntermediateLayer.inMeter;
            pos = psfObj.position.inMeter;
            z = pos(3);
            focalLength = psfObj.objectiveFocalLength.inMeter;
            mu = 1e-12;

            %% Pre-Calculations

            if length(RI)==1
                RI=[RI, RI, RI];
            end

            % Coordinates in the objective pupil
            par.nGrid = obj.nGrid;
            par.spacing = obj.unitKSpace;
            pupilMask = Mask(par);
            [~, maskAngle] = pupilMask.getPolarCoordinates;
            PHI3 = fliplr(maskAngle) - pi;

            % Wavenumbers (magnitude of k-vectors) in the different media, k = k0 * RI
            k0 = 2 * pi / psfObj.wavelength.inMeter; % wavenumber in vacuum
            k1 = k0 * RI(1); % wavenumber in media 1 (typically water), Eq. (14) of paper
            k2 = k0 * RI(2); % wavenumber in media 2 (intermediate layer), Eq. (14) of paper
            k3 = k0 * RI(3); % wavenumber in media 3 (immersion medium), Eq. (14) of paper

            % Angles in different media
            Kr = pupilMask.radius;
            THETA1 = acos( sqrt( 1 - (RI(3)/RI(1) * Kr/k3).^2 ) ) .* pupilMask.values; % angle in medium 1; Eq. (4) of paper
            THETA2 = acos( sqrt( 1 - (RI(3)/RI(2) * Kr/k3).^2 ) ) .* pupilMask.values; % angle in medium 2; Eq. (4) of paper
            THETA3 = asin( Kr/k3 ) .* pupilMask.values; % angle in medium 3; maximum theta3 from Eq. (19) of paper

            %% Calculations according to paper of Axelrod, 2012

            % Cosines of angles
            CTHETA1 = cos(THETA1);
            CTHETA2 = cos(THETA2);
            CTHETA3 = cos(THETA3);

            % Fresnel-coefficients
            % Eq. (3) of paper
            tp12 = 2*RI(1)*CTHETA1./(RI(1)*CTHETA2+RI(2)*CTHETA1);
            tp23 = 2*RI(2)*CTHETA2./(RI(2)*CTHETA3+RI(3)*CTHETA2);

            ts12 = 2*RI(1)*CTHETA1./(RI(1)*CTHETA1+RI(2)*CTHETA2);
            ts23 = 2*RI(2)*CTHETA2./(RI(2)*CTHETA2+RI(3)*CTHETA3);

            rp12 = (RI(2)*CTHETA1-RI(1)*CTHETA2)./(RI(1)*CTHETA2+RI(2)*CTHETA1);
            rp23 = (RI(3)*CTHETA2-RI(2)*CTHETA3)./(RI(2)*CTHETA3+RI(3)*CTHETA2);

            rs12 = (RI(1)*CTHETA1-RI(2)*CTHETA2)./(RI(1)*CTHETA1+RI(2)*CTHETA2);
            rs23 = (RI(2)*CTHETA2-RI(3)*CTHETA3)./(RI(2)*CTHETA2+RI(3)*CTHETA3);

            % Fresnel coefficients for three-layer system
            % Eq. (12) of paper
            tp = tp12 .* tp23 .* exp(1i*k2*hIntermediate*CTHETA2) ./ (1 + rp12 .* rp23 .* exp(2i*k2*hIntermediate*CTHETA2));
            ts = ts12 .* ts23 .* exp(1i*k2*hIntermediate*CTHETA2) ./ (1 + rs12 .* rs23 .* exp(2i*k2*hIntermediate*CTHETA2));

            % Dipole projections onto directions p, s and z
            % Eq. (13) of paper
            mu_p = mu * sin(dipole.inclination) .* cos(dipole.azimuth - PHI3);
            mu_s = mu * sin(dipole.inclination) .* sin(dipole.azimuth - PHI3);
            mu_z = mu * cos(dipole.inclination);

            % Prefactor C (the only constant where f plays a role)
            % Eq. (11) of paper
            C = ( k3^2 * exp(1i*k3*focalLength) .* CTHETA3) / (focalLength * RI(1)) ...
                .* exp(-1i*k3*hIntermediate*CTHETA3) ...
                .* exp(1i*k1.*CTHETA1.*z);

            % Electric field components in layer 3 (pre-objective zone), along the p, s and z-axes
            % (see paper for axes definitions)
            % Eq. (10) of paper
            % See Erratum of paper for correct equations (https://onlinelibrary.wiley.com/doi/10.1111/jmi.12173)
            % Previously, the factor n1 in the denominator of the second term in the parenthesis was missing in the first and third line
            E3p = C .* tp .* CTHETA3 .* (mu_p./RI(3) + mu_z.*sin(THETA3)./(RI(1)*CTHETA1));
            E3s = C .* ts .* (mu_s./(RI(3)./CTHETA1));
            E3z = C .* tp .* sin(THETA3) .* (mu_p./RI(3) + mu_z.*sin(THETA3)./(RI(1)*CTHETA1));

            % Influence of objective, rotation of rays by their angle theta3 such that they are all parallel to the optical axis
            % Eq. (15) and (16) of paper
            apodizationFactor = 1 ./ sqrt(CTHETA3) .* pupilMask.values; % apodization of objective lens
            E_BFP_p = (E3p.*CTHETA3 + E3z.*sin(THETA3)) .* apodizationFactor;
            E_BFP_s = E3s .* apodizationFactor;  % s-polarization remains unchanged by this rotation

            % Coordinate transformation into x-and y-polarization yields fields in the back focal plane of objective
            % Eq. (18) of paper
            E_BFP.x = cos(PHI3).*E_BFP_p - sin(PHI3).*E_BFP_s;
            E_BFP.y = sin(PHI3).*E_BFP_p + cos(PHI3).*E_BFP_s;

        end

        function plot(obj)
            %% Electric field
            % Electric field x
            subplot(2,2,1)
            imagesc(real(obj.electricField.x))
            axis equal; axis tight
            title('E_x real part')

            subplot(2,2,2)
            imagesc(imag(obj.electricField.x))
            axis equal; axis tight
            title('E_x imaginary part')

            % Electric field y
            subplot(2,2,3)
            imagesc(real(obj.electricField.y))
            axis equal; axis tight
            title('E_y real part')

            subplot(2,2,4)
            imagesc(imag(obj.electricField.y))
            axis equal; axis tight
            title('E_y imaginary part')

            %% Intensity
            figure
            subplot(1,3,1)
            imagesc(abs(obj.electricField.x).^2)
            axis equal; axis tight
            title('Intensity x-pol.')

            subplot(1,3,2)
            imagesc(abs(obj.electricField.y).^2)
            axis equal; axis tight
            title('Intensity y-pol.')

            subplot(1,3,3)
            imagesc(abs(obj.electricField.x).^2 + abs(obj.electricField.y).^2)
            axis equal; axis tight
            title('Total intensity')
        end
    end
end
