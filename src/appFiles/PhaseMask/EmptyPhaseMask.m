classdef EmptyPhaseMask < PhaseMask

    methods
        function obj = EmptyPhaseMask(gridSize)
            obj.mask = zeros(gridSize);
        end
    end

    methods (Access=public)
        % Overwrite apply method:
        function shiftedElectricField = apply(obj, electricField)
            shiftedElectricField = electricField;
        end
    end
end