classdef IsotropicPSF < PSF

    methods
        %% Functions for calculation of PSF
        function obj = createImage(obj)
            [obj.Defocus, obj.pupilMask, obj.chirpZTransform, obj.phaseMaskObj, ...
                obj.zernikeAberrationsObj, obj.attenuationMaskObj, obj.transmissionMaskObj] = obj.setup();
            
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
            midSlice = ceil(size(psf,3)/2);
            totalIntensity = sum(psf(:,:,midSlice), 'all');
            psf = psf ./ totalIntensity * obj.nPhotons;

            % Ix
            Ix = Ixx  + Ixy + Ixz;
            Ix = Ix .* repmat(obj.pixelSensitivityMask, obj.nPixels);
            Ix = reduceOversampling(obj, Ix);
            Ix = Ix ./ totalIntensity * obj.nPhotons;

            % Iy
            Iy = Iyx  + Iyy + Iyz;
            Iy = Iy .* repmat(obj.pixelSensitivityMask, obj.nPixels);
            Iy = reduceOversampling(obj, Iy);
            Iy = Iy ./ totalIntensity * obj.nPhotons;

            psf = adjustExcitation(obj, psf);
            psf = applyShotNoise(obj, psf);
            psf = addBackgroundNoise(obj, psf);
            I_x = applyShotNoise(obj, Ix);
            I_y = applyShotNoise(obj, Iy);
            obj.image = psf;
            obj.Ix = I_x;
            obj.Iy = I_y; 
        end

        function [Ix,Iy] = getDipoleImage(obj,dipole)
            obj.dipole = dipole;
            bfp = BackFocalPlane(obj);
            obj.backFocalPlane = bfp;

            % Apply phase mask
            obj.fieldBFP.x = obj.phaseMaskObj.apply(bfp.electricField.x);
            obj.fieldBFP.y = obj.phaseMaskObj.apply(bfp.electricField.y);

            % Apply Zernike aberrations mask
            obj.fieldBFP.x = obj.zernikeAberrationsObj.apply(obj.fieldBFP.x);
            obj.fieldBFP.y = obj.zernikeAberrationsObj.apply(obj.fieldBFP.y);

            % Apply attenuation mask
            obj.fieldBFP.x = obj.attenuationMaskObj.apply(obj.fieldBFP.x);
            obj.fieldBFP.y = obj.attenuationMaskObj.apply(obj.fieldBFP.y);

            % Apply transmission mask
            obj.fieldBFP.x = obj.transmissionMaskObj.apply(obj.fieldBFP.x);
            obj.fieldBFP.y = obj.transmissionMaskObj.apply(obj.fieldBFP.y);

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

