classdef Mask

    properties
        values double
        nGrid = 129
        spacing = 1 % spacing between two neighbouring entries
        mode = 'FFT' % 'exact' or 'FFT'
                     % 'exact':
                     % zero point is exactly in the middle,
                     % i.e. for even grid sizes between the two middle pixels
                     % 'FFT':
                     % zero point is the middle pixel for odd grid sizes
                     % and the pixel to the lower right of the exact centre for even grid sizes
                     % (this is also the pixel which MATLAB takes as zero)
        radius
    end

    methods
        function obj = Mask(par)
            if nargin > 0
                obj = setInputParameters('Mask', obj, par);
            end

            % Calculate mask values
            N = [obj.nGrid, obj.nGrid];
            s = [obj.spacing, obj.spacing];
            if strcmp(obj.mode,'exact')==1
                x = linspace( -(N(1)-1)/2, (N(1)-1)/2, N(1) ) * s(1);
                y = linspace( -(N(2)-1)/2, (N(2)-1)/2, N(2) ) * s(2);
            elseif strcmp(obj.mode,'FFT')==1
                x = ( floor( -N(1)/2 + 0.5) : floor( N(1)/2 - 0.5) ) * s(1);
                y = ( floor( -N(2)/2 + 0.5) : floor( N(2)/2 - 0.5) ) * s(2);
            end

            [X, Y] = meshgrid(y,x);
            obj.radius = sqrt(X.^2+Y.^2);
            obj.values = ( obj.radius.^2 <= ((min(N/2.*s))).^2 );
        end

        function [normalizedRadius, angle] = getPolarCoordinates(obj)
            x = linspace(-1,1,obj.nGrid);
            [X,Y] = meshgrid(x,x);
            [angle,normalizedRadius] = cart2pol(X,Y);
            angle = fliplr(angle + pi);

            %idxNan = (obj.values==0);
            %angle(idxNan) = 0; %NaN;
            %radius(idxNan) = 0; %NaN;
        end

        function plot(obj)
            imagesc(obj.values)
            axis equal; axis tight;
        end
    end
end