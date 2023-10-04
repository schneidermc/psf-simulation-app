classdef PyramidPhaseMask < PhaseMask

    properties
        nFacets {mustBeInteger, mustBeGreaterThan(nFacets,1)} = 4
    end
    
    methods
        function obj = PyramidPhaseMask(gridSize, nFacets, maxShift)
            if nargin < 3
                maxShift = 2*pi;
            end
            
            obj.nFacets = nFacets;
            
            sectorAngles = [0, 2*pi/obj.nFacets] - 2*pi/obj.nFacets/2 + pi/2;
            combinedObj = ConstantPhaseMask(gridSize, 0);

            for k = 1:obj.nFacets
                thisSectorAngles = sectorAngles;
                thisSector = TiltPhaseMask(gridSize, maxShift).selectCircleSector(thisSectorAngles);
                thisSector = thisSector.rotate(mod(-pi/2 + (k-1)*2*pi/obj.nFacets,2*pi));
                combinedObj = (combinedObj + thisSector);
            end

            obj.mask = combinedObj.mask;
        end
    end
end