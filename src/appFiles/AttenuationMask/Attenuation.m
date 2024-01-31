classdef Attenuation
    
    properties
        mask (:,:) double {mustBeInRange(mask,0,1)}
    end

    properties (Dependent, Hidden)
        gridSize
        polarCoordinates
    end


    methods
        %% Constructor
        function obj = Attenuation(mask)
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
            [theta,rho] = Attenuation.createGrid(obj.gridSize);
            polarCoord.angle = theta;
            polarCoord.radius = rho;
        end
    end

    methods (Access=public)
        %% Functions
        function attenuatedElectricField = apply(obj, electricField)
            assert(isequal(size(obj.mask),size(electricField)))
            attenuatedElectricField = sqrt(1-obj.mask) .* electricField;
            % Note: Full attenuation (=1) ~= no electric field (.*0)
            % Attenuation coefficient is for intensity;
            % take sqrt to apply attenuation directly to electric field
        end

        %% Plot
        function plot(obj, axesHandle)
            if nargin < 2
                axesHandle = gca;
            end
            imagesc(axesHandle, obj.mask)
            axis(axesHandle, 'equal');
            axis(axesHandle, 'tight');
            cb = colorbar(axesHandle);
            colormap(axesHandle, flipud(gray))
            caxis(axesHandle, [0 1])
            cb.Label.String = 'Attenuation';
            cb.Label.FontSize = 12;
            cb.Ticks = (0:0.25:1);
            cb.TickLabels = {'0','','0.5','','1'};
            set(axesHandle,'visible','off')
        end
    end

    methods (Access=public)
        %% Adapt attenuation mask
        function obj = cutInnerRing(obj, relativeRadius)
            mustBeInRange(relativeRadius,0,1)
            obj.mask(obj.polarCoordinates.radius < relativeRadius) = 1;
        end

        function obj = cutAperture(obj, relativeRadius)
            if nargin < 2
                relativeRadius = 1;
            end
            mustBeInRange(relativeRadius,0,1)
            obj.mask(obj.polarCoordinates.radius > relativeRadius) = 1;
        end

        function obj = selectCircleSector(obj, sectorAngles, attenuation)
            mustBeInFullRadialRange(sectorAngles)
            if numel(sectorAngles)==1
                sectorAngles = [0,sectorAngles];
            end
            if nargin < 3
                attenuation = 1;
            end
            isInSegment = ((sectorAngles(1)<obj.polarCoordinates.angle)) & (obj.polarCoordinates.angle<sectorAngles(2));
            obj.mask(isInSegment) = attenuation;
            obj = obj.cutAperture();
        end

        function obj = rotate(obj,rotationAngle)
            mustBeInFullRadialRange(rotationAngle)
            if rotationAngle ~= 0
                obj.mask = imrotate(obj.mask, -rotationAngle*180/pi, 'bilinear', 'crop');
                % Fill cropped area with ones
                obj.mask(obj.polarCoordinates.radius > 1) = 1;
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