classdef Transmission

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
        function obj = Transmission(mask)
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
        function attenuatedElectricField = apply(obj, electricField)
            assert(isequal(size(obj.mask),size(electricField)))
            attenuatedElectricField = obj.mask .* electricField; % apply amplitude modulation
        end

        %% Plot
        function plot(obj, axesHandle)
            if nargin < 2
                axesHandle = gca;
            end
            h = imagesc(axesHandle, obj.mask);
            set(h, 'AlphaData', ~isnan(obj.pupil))
            axis(axesHandle, 'equal');
            axis(axesHandle, 'tight');
            cb = colorbar(axesHandle);
            colormap(axesHandle, gray)
            caxis(axesHandle, [0 1])
            cb.Label.String = 'Transmission';
            cb.Label.FontSize = 16;
            cb.Ticks = (0:0.25:1);
            cb.TickLabels = {'0','','0.5','','1'};
            set(axesHandle,'visible','off')
        end
    end

    methods (Access=public)
        %% Adapt transmission mask
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
    end
    
    methods (Static)
        function [theta,rho] = createGrid(gridSize)
            x = linspace(-1, 1, gridSize);
            [theta, rho] = cart2pol(x', x);
            theta = flipud(mod(theta,2*pi));
        end
    end
end
