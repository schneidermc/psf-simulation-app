classdef Vortex < PhaseMask

    methods
        function obj = Vortex(gridSize)
            [theta, ~] = PhaseMask.createGrid(gridSize);
            obj.mask = theta;
            obj = obj.cutAperture();
        end
    end
end