classdef ConstantAttenuationMask < Attenuation
    
    methods
        function obj = ConstantAttenuationMask(gridSize, attenuation)
            obj.mask = attenuation*ones(gridSize);
            obj = obj.cutAperture();
        end
    end
end