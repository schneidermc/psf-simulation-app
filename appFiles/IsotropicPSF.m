classdef IsotropicPSF < PSF

    properties
        % Fluorophore
        rotation (1,1) logical = 0
    end

    methods
        %% Functions for calculation of PSF
        function obj = createImage(obj)
            [obj.Defocus, obj.pupilMask, obj.chirpZTransform, obj.phaseMaskObj] = obj.setup();

            % xyz-dipoles
            xDipole = Dipole(pi/2,0);
            [Ixx,Iyx] = getDipoleImage(obj,xDipole);
            yDipole = Dipole(pi/2,pi/2);
            [Ixy,Iyy] = getDipoleImage(obj,yDipole);
            zDipole = Dipole(0,0);
            [Ixz,Iyz] = getDipoleImage(obj,zDipole);

            psf = Ixx + Iyx + Ixy + Iyy + Ixz + Iyz;

            % Account for pixel sensitivity
            psf = psf .* repmat(obj.pixelSensitivityMask, obj.nPixels);

            % Sum over block matrices (to return to desired pixelsize)
            psf = reduceOversampling(obj, psf);

            % Normalization
            totalIntensity = sum(sum(psf));
            psf = psf ./ totalIntensity * obj.nPhotons;

            psf = adjustExcitation(obj, psf);
            psf = applyShotNoise(obj, psf);
            psf = addBackgroundNoise(obj, psf);
            obj.image = psf;
        end

        function [Ix,Iy] = getDipoleImage(obj,dipole)
            obj.dipole = dipole;
            bfp = BackFocalPlane(obj);
            obj.backFocalPlane = bfp;

            % Apply phase mask
            obj.fieldBFP.x = obj.phaseMaskObj.apply(bfp.electricField.x);
            obj.fieldBFP.y = obj.phaseMaskObj.apply(bfp.electricField.y);

            n = obj.nPixels*obj.oversampling;
            Ix = zeros(n,n);
            Iy = zeros(n,n);
            for k=1:size(obj.stageDrift.motion,1)
                % Apply aberrations
                aberrations = getAberrations(obj,k);
                aberratedFieldBFP = applyAberrations(obj, aberrations);

                % Get intensities from BFP field
                [Ix_tmp,Iy_tmp] = getIntensitiesFromBFP(obj, aberratedFieldBFP);
                Ix = Ix + Ix_tmp;
                Iy = Iy + Iy_tmp;
            end
        end


        function [I_x, I_y] = getIntensitiesFromBFP(obj, fieldBFP)
            I_x = abs( obj.chirpZTransform.apply(obj, fieldBFP.x) ).^2; % intensity = |E_imagePlane|²
            I_y = abs( obj.chirpZTransform.apply(obj, fieldBFP.y) ).^2; % intensity = |E_imagePlane|²
        end

        %% Plot function
        function plot(obj,setLimitsFullRange)
            if nargin==2 && setLimitsFullRange
                imagesc(obj.image(:,:),[0 max(obj.image(:))]);
            else
                imagesc(flip(obj.image(:,:)));
            end
            colorbar;
            axis equal; axis tight; axis xy;
            title(['Pixelsize = ', num2str(obj.pixelSize.inNanometer), 'nm. ', ...
                'Rotating Dipole']);
        end
    end
end

