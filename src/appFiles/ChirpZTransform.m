classdef ChirpZTransform
    % 2D chirp-z transform
    %
    % See the following publication:
    % Raoqiong Tong and Robert W. Cox.
    % "Rotation of NMR images using the 2D chirp-z transform"
    % Magnetic resonance in medicine 41.2 (1999): 253-256.
    
    properties
        kernel
        fourierKernel
        nPixelsPadded
    end

    methods
        function obj = ChirpZTransform(psfObj)
            % Get dimensions of image and k-space
            nPixelsImage = psfObj.nPixels .* psfObj.oversampling; % size of field in pixel
            nPixelsBFP = psfObj.nDiscretizationBFP; % size of input field (k-space)
            obj.nPixelsPadded = nPixelsImage + nPixelsBFP - 1; % padded size required for convolution is Nx+Nk-1

            % Create grid
            x = ( floor(-obj.nPixelsPadded(1)/2 + 0.5) : floor(obj.nPixelsPadded(1)/2 - 0.5) )';
            y = ( floor(-obj.nPixelsPadded(1)/2 + 0.5) : floor(obj.nPixelsPadded(1)/2 - 0.5) );

            % Create kernel (quadratic convolution phase kernel)
            alpha = psfObj.unitObjectSpace * psfObj.unitKSpace / (2*pi);
            obj.kernel = exp(-1i*alpha*pi*x.^2) * exp(-1i*alpha*pi*y.^2); % = exp(-1i*alpha*pi*(x.^2+y.^2)), but slightly faster
            
            % Fourier transform of kernel
            obj.fourierKernel = fft2( ifftshift(obj.kernel) );
        end

        function E_out = apply(obj,psfObj,E_in)
            % Complex conjugate of kernel
            conjKernel = conj(obj.kernel);
            
            % Electric field times complex conjugate of kernel
            if ismatrix(E_in)
                f = embedArray2D(E_in,obj.nPixelsPadded,0) .* conjKernel;
            else
                f = embedArray3D(E_in,obj.nPixelsPadded,0) .* conjKernel;
            end

            % Fourier transform of kernel
            F = fft2( ifftshift(f) );

            % Convolving f with kernel (inverse fft of multiplication)
            convolution = fftshift( ifft2( F .* obj.fourierKernel ) );

            % Final multiplication by factor
            E_out = obj.nPixelsPadded * conjKernel .* convolution;

            % Crop to image size
            if ismatrix(E_in)
                E_out = cropArray2D(E_out, psfObj.nPixels .* psfObj.oversampling);
            else
                E_out = cropArray3D(E_out, psfObj.nPixels .* psfObj.oversampling);
            end
        end
    end
end