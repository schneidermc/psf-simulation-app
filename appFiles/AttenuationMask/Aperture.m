classdef Aperture < Attenuation
    
    methods
        function obj = Aperture(gridSize)
            sector = ConstantAttenuationMask(gridSize,0);
            obj.mask = sector.mask;
        end
    end
end