classdef TiltPhaseMask < PhaseMask
    
    methods
        function obj = TiltPhaseMask(gridSize, maxShift)
            if nargin < 2
                maxShift = 2*pi;
            end

            x = linspace(-maxShift,maxShift,gridSize);
            [X,~] = meshgrid(x,x);
            X = abs(X);
            
            obj.mask = X;
            obj = obj.cutAperture();
        end
    end
end