classdef PhaseMask

    properties
        mask (:,:) double
    end

    properties (Dependent, Hidden)
        gridSize
        polarCoordinates
        pupil
    end


    methods
        %% Constructor
        function obj = PhaseMask(mask)
            if nargin > 0
                obj.mask = mask;
            end
        end

        %% Get methods
        function mask = getMask(obj)
            mask = obj.mask;
        end
        function gridSize = get.gridSize(obj)
            gridSize = size(obj.mask,1);
        end
        function polarCoord = get.polarCoordinates(obj)
            [theta,rho] = PhaseMask.createGrid(obj.gridSize);
            polarCoord.angle = theta;
            polarCoord.radius = rho;
        end
        function pupil = get.pupil(obj)
            parMask.nGrid = obj.gridSize;
            pupilMask = Mask(parMask);
            values = pupilMask.values;
            values(values==0) = NaN;
            pupil = values;
        end
    end

    methods (Access=public)
        %% Functions
        function shiftedElectricField = apply(obj, electricField)
            assert(isequal(size(obj.mask),size(electricField)))
            shiftedElectricField = electricField .* exp(1i.*obj.mask);
        end

        %% Plot
        function plot(obj, axesHandle)
            if nargin < 2
                axesHandle = gca;
            end
            h = imagesc(axesHandle, mod(obj.mask,2*pi));
            set(h, 'AlphaData', ~isnan(obj.pupil))
            axis(axesHandle, 'equal');
            axis(axesHandle, 'tight');
            cb = colorbar(axesHandle);
            caxis(axesHandle, [0 2*pi])
            cb.Label.String = 'Phase shift';
            cb.Label.FontSize = 16;
            cb.Ticks = (0:0.5:2)*pi;
            cb.TickLabels = {'0','','\pi','','2\pi'};
            set(axesHandle,'visible','off')
        end
    end

    methods (Access=public)
        %% Adapt phase mask
        function obj = cutInnerRing(obj, relativeRadius)
            mustBeInRange(relativeRadius,0,1)
            obj.mask(obj.polarCoordinates.radius < relativeRadius) = 0;
        end

        function obj = cutAperture(obj, relativeRadius)
            if nargin < 2
                relativeRadius = 1;
            end
            mustBeInRange(relativeRadius,0,1)
            if relativeRadius == 1
                values = obj.pupil;
                values(isnan(values)) = 0;
                obj.mask = obj.mask .* values;
            else
                obj.mask(obj.polarCoordinates.radius > relativeRadius) = 0;
            end
        end

        function obj = selectCircleSector(obj, sectorAngles)
            mustBeInFullRadialRange(sectorAngles)
            if numel(sectorAngles)==1
                sectorAngles = [0,sectorAngles];
            end
            isInSegment = ((sectorAngles(1)<obj.polarCoordinates.angle)) & (obj.polarCoordinates.angle<sectorAngles(2));
            obj.mask(~isInSegment) = 0;
        end

        function obj = rotate(obj,rotationAngle)
            mustBeInFullRadialRange(rotationAngle)
            if rotationAngle ~= 0
                obj.mask = imrotate(obj.mask, -rotationAngle*180/pi, 'bilinear', 'crop');
            end
        end

        %% Overload operators
        function combinedPhaseMask = plus(obj1,obj2)
            combinedPhaseMask = PhaseMask(obj1.mask + obj2.mask);
        end

        function combinedPhaseMask = times(obj1,obj2)
            combinedPhaseMask = PhaseMask(obj1.mask .* obj2.mask);
        end
    end

    methods (Static)
        function [theta,rho] = createGrid(gridSize)
            x = linspace(-1, 1, gridSize);
            [theta, rho] = cart2pol(x', x);
            theta = flipud(mod(theta,2*pi));
        end
    end
end