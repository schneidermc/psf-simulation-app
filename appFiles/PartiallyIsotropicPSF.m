classdef PartiallyIsotropicPSF < IsotropicPSF
    
    methods
        function obj = PartiallyIsotropicPSF(par)
            obj = obj@IsotropicPSF(par);
            obj.rotationalConstraint = par.rotationalConstraint;
        end 

        %% Functions for calculation of PSF
        function obj = createImage(obj)
            [obj.Defocus, obj.pupilMask, obj.chirpZTransform, obj.phaseMaskObj, ...
                obj.zernikeAberrationsObj, obj.attenuationMaskObj, obj.transmissionMaskObj] = obj.setup();
            
            % Fixed dipole
            [Ix_fixed,Iy_fixed] = getDipoleImage(obj,obj.dipole);
            % xyz-dipoles
            xDipole = Dipole(pi/2,0);
            [Ixx,Iyx] = getDipoleImage(obj,xDipole);
            yDipole = Dipole(pi/2,pi/2);
            [Ixy,Iyy] = getDipoleImage(obj,yDipole);
            zDipole = Dipole(0,0);
            [Ixz,Iyz] = getDipoleImage(obj,zDipole);

            psf = (1 - obj.rotationalConstraint)/3 * (Ixx + Iyx + Ixy + Iyy + Ixz + Iyz) + obj.rotationalConstraint/3 * (Ix_fixed + Iy_fixed);

            % Account for pixel sensitivity
            psf = psf .* repmat(obj.pixelSensitivityMask, obj.nPixels);

            % Sum over block matrices (to return to desired pixelsize)
            psf = reduceOversampling(obj, psf);

            % Normalization
            totalIntensity = sum(sum(psf));
            psf = psf ./ totalIntensity * obj.nPhotons;

            % Ix
            Ix = (1 - obj.rotationalConstraint)/3 * (Ixx  + Ixy + Ixz) + obj.rotationalConstraint/3 * Ix_fixed;
            Ix = Ix .* repmat(obj.pixelSensitivityMask, obj.nPixels);
            Ix = reduceOversampling(obj, Ix);
            Ix = Ix ./ totalIntensity * obj.nPhotons;

            % Iy
            Iy = (1 - obj.rotationalConstraint)/3 * (Iyx  + Iyy + Iyz) + obj.rotationalConstraint/3 * Iy_fixed;
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
                'Partially Rotating Dipole']);
        end
    end
end

