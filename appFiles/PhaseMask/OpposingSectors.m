classdef OpposingSectors < PhaseMask
    
    methods
        function obj = OpposingSectors(gridSize, sectorAngles, constantShift)
            if nargin < 3
                constantShift = pi;
            end
            sector1 = ConstantPhaseMask(gridSize,constantShift).selectCircleSector(sectorAngles);
            sector2 = ConstantPhaseMask(gridSize,constantShift).selectCircleSector(sectorAngles+pi);
            combinedObj = (sector1 + sector2);
            obj.mask = combinedObj.mask;
        end
    end
end