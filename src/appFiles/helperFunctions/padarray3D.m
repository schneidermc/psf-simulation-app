function arrayOut = padarray3D(arrayIn, targetSize, padVal)
    % Embeds a 3D array inside a bigger canvas with size targetSize.
    % The frame around is filled with padval.
    % targetSize ... size of arrayOut
    % padVal ... padding value

    arguments
        arrayIn (:,:,:)
        targetSize (1,:) {mustBeInteger, mustBePositive} = size(arrayIn)
        padVal (1,1) = 0
    end

    if length(targetSize)==1
        targetSize = [targetSize, targetSize];
    end

    sizeIn = size(arrayIn);
    sizeDiff = targetSize(1:2) - sizeIn(1:2); % [nRows, nCols] to be added
    
    if any(sizeDiff<0)
        error('Target size for embedding must be greater than size of input array!')
    else
        padSize = [ceil(sizeDiff(1)/2), floor(sizeDiff(1)/2)];
        arrayOut = padarray(arrayIn, padSize, padVal, 'both');
    end
end
