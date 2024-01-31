classdef PixelSensitivity
    
    properties
        mask
    end
    
    enumeration
        centerLarge (PixelSensitivity.centerMask(9,1))
        centerSmall (PixelSensitivity.centerMask(9,3))
        centerOnly (PixelSensitivity.centerMask(3,1))
        none (0)
    end
    
    methods
        function obj = PixelSensitivity(mask)
            obj.mask = mask;
        end
    end
    
    methods (Static)
        function mask = uniform(n)
            mask = ones(n);
        end
        
        function mask = centerMask(width, border)
            mask = ones(width);
            mask(:,1:border) = 0;
            mask(:,end-border+1:end) = 0;
            mask(1:border,:) = 0;
            mask(end-border+1:end,:) = 0;
        end
    end
end

