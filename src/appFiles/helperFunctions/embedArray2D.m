function arrayOut = embedArray2D(arrayIn, targetSize, padVal)
    % Embeds a 2D array inside a bigger canvas with size targetSize.
    % The frame around is filled with padval.
    % targetSize ... size of arrayOut
    % padVal ... padding value

    arguments
        arrayIn (:,:)
        targetSize (1,:) {mustBeInteger, mustBePositive} = size(arrayIn)
        padVal (1,1) = 0
    end

    if length(targetSize)==1
        targetSize = [targetSize, targetSize];
    end

    sizeIn = size(arrayIn);
    sizeDiff = targetSize - sizeIn; % [nRows, nCols] to be added
    
    if any(sizeDiff<0)
        error('Target size for embedding must be greater than size of input array!')
    else
        arrayOutTmp = simplePadAboveBelow(arrayIn, [ceil(sizeDiff(1)/2), floor(sizeDiff(1)/2)], padVal);
        arrayOut = simplePadLeftRight(arrayOutTmp, [ceil(sizeDiff(2)/2), floor(sizeDiff(2)/2)], padVal);
    end
end

function arrayOut = simplePadAboveBelow(arrayIn, nAboveBelow, padval)
    width = size(arrayIn, 2);
    arrayOut = [padval.*ones(nAboveBelow(1),width); arrayIn; padval.*ones(nAboveBelow(2),width)];
end

function arrayOut = simplePadLeftRight(arrayIn, nLeftRight, padval)
    height = size(arrayIn, 1);
    arrayOut = [padval.*ones(height, nLeftRight(1)), arrayIn, padval.*ones(height, nLeftRight(2))];
end