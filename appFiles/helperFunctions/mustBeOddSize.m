function mustBeOddSize(a)
    n = size(a,1);
    if ~(rem(n,2) == 1)
        eidType = 'mustBeOdd:notOdd';
        msgType = 'Input must be an odd number.';
        throwAsCaller(MException(eidType,msgType))
    end
end