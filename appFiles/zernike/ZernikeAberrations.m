classdef ZernikeAberrations < PhaseMask

    properties
        zernikeNollIndices (1,:) = []
        zernikeCoefficients (1,:) = []
    end

    methods
        %% Constructor
        function obj = ZernikeAberrations(zernikeNollIndices, zernikeCoefficients, nGrid)
            if nargin < 3
                nGrid = 129;
            end
            if nargin > 0
                parMask.nGrid = nGrid;
                pupilMask = Mask(parMask);
                Z = ZernikePolynomials.getInstance(pupilMask);
                aberrations = Z.getAberration(zernikeNollIndices,zernikeCoefficients);
                obj.mask = 2*pi*aberrations;
            end
        end
        
        %% Plot
        function plot(obj, axesHandle)
            if nargin < 2
                axesHandle = gca;
            end
            h = imagesc(axesHandle, (mod(obj.mask+pi,2*pi)-pi));
            set(h, 'AlphaData', ~isnan(obj.pupil))
            axis(axesHandle, 'equal');
            axis(axesHandle, 'tight');
            cb = colorbar(axesHandle);
            caxis(axesHandle, [-pi pi])
            cb.Label.String = 'Phase shift';
            cb.Label.FontSize = 16;
            cb.Ticks = (-1:0.5:1)*pi;
            cb.TickLabels = {'-\pi','','0','','\pi'};
            set(axesHandle,'visible','off')
        end
    end
end