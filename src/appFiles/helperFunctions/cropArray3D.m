function arrayOut = cropArray3D(arrayIn, targetSize)
    % Crops a 3D array to size targetSize.
    % Cropping is done with respect to center.
    % targetSize ... size of arrayOut

    arguments
        arrayIn (:,:,:)
        targetSize (1,:) {mustBeInteger, mustBePositive} = size(arrayIn)
    end

    if length(targetSize)==1
        targetSize = [targetSize, targetSize];
    end

    if length(targetSize)==2
        targetSize = [targetSize, size(arrayIn, 3)];
    end
    
    sizeIn = size(arrayIn);
    sizeDiff = sizeIn(1:2) - targetSize(1:2); % [nRows, nCols] to be removed
    
    if any(sizeDiff<0)
        error('Target size for cropping must be smaller than size of input array!')
    else
        indKeepRow = (floor(sizeDiff(1)/2) + 1) : (floor(sizeDiff(1)/2) + targetSize(1));
        indKeepCol = (floor(sizeDiff(2)/2) + 1) : (floor(sizeDiff(2)/2) + targetSize(2));
        arrayOut = arrayIn(indKeepRow,indKeepCol,:);
    end
end
