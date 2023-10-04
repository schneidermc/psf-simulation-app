function arrayOut = embedArray2DFill(arrayIn, targetSize, padVal)
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
    
    % Create whole array of padVal with right size
    arrayOut = padVal .* ones(targetSize);
    
    % Fill in input array
    indRow = (floor(sizeDiff(1)/2) + 1) : (targetSize(1) - floor(sizeDiff(1)/2));
    indCol = (floor(sizeDiff(2)/2) + 1) : (targetSize(2) - floor(sizeDiff(2)/2));
    arrayOut(indRow,indCol) = arrayIn;
end