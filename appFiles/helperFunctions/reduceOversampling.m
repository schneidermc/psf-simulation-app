function [psfOut, IxOut, IyOut] = reduceOversampling(obj, psf, Ix, Iy)
    if ismatrix(psf)
        psfOut = squeeze(sum(sum(reshape(psf,obj.oversampling,obj.nPixels,obj.oversampling,obj.nPixels),1),3));
        if nargin > 2
            IxOut = squeeze(sum(sum(reshape(Ix,obj.oversampling,obj.nPixels,obj.oversampling,obj.nPixels),1),3));
            IyOut = squeeze(sum(sum(reshape(Iy,obj.oversampling,obj.nPixels,obj.oversampling,obj.nPixels),1),3));
        else
            IxOut = [];
            IyOut = [];
        end
    elseif ndims(psf)
        psfOut = squeeze(sum(sum(reshape(psf,obj.oversampling,obj.nPixels,obj.oversampling,obj.nPixels,size(psf,3)),1),3));
        if nargin > 2
            IxOut = squeeze(sum(sum(reshape(Ix,obj.oversampling,obj.nPixels,obj.oversampling,obj.nPixels,size(Ix,3)),1),3));
            IyOut = squeeze(sum(sum(reshape(Iy,obj.oversampling,obj.nPixels,obj.oversampling,obj.nPixels,size(Iy,3)),1),3));
        else
            IxOut = [];
            IyOut = [];
        end
    else
        error('Input data must be 2D or 3D')
    end
end