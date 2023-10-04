classdef Astigmatism < PhaseMask

    methods
        function obj = Astigmatism(gridSize, astigmatismCoefficient)
            par.nGrid = gridSize;
            pupil = Mask(par);
            Z = ZernikePolynomials.getInstance(pupil);
            aberration = Z.getAberration(6,astigmatismCoefficient);
            obj.mask = aberration;
        end

        function plot(obj, axesHandle)
            if nargin < 2
                axesHandle = gca;
            end
            imagesc(axesHandle, obj.mask)
            axis(axesHandle, 'equal');
            axis(axesHandle, 'tight');
            cb = colorbar(axesHandle);
            caxis(axesHandle, [min(obj.mask(:)) max(obj.mask(:))])
            cb.Label.String = 'Phase shift';
            cb.Label.FontSize = 12;
            set(axesHandle,'visible','off')
        end
    end
end