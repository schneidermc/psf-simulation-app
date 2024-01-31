classdef Sector < PhaseMask
    
    methods
        function obj = Sector(gridSize, sectorAngles, constantShift)
            if nargin < 3
                constantShift = pi;
            end
            sector = ConstantPhaseMask(gridSize,constantShift).selectCircleSector(sectorAngles);
            obj.mask = sector.mask;
        end
    end
end