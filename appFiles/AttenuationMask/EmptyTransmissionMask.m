classdef EmptyTransmissionMask < Transmission

    methods
        function obj = EmptyTransmissionMask(gridSize)
            obj.mask = ones(gridSize);
            obj = obj.cutAperture();
        end
    end

    methods (Access=public)
        % Overwrite apply method:
        function shiftedElectricField = apply(obj, electricField)
            shiftedElectricField = electricField;
        end
    end
end
