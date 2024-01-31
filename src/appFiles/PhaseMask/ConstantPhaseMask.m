classdef ConstantPhaseMask < PhaseMask
    
    methods
        function obj = ConstantPhaseMask(gridSize, constantShift)
            obj.mask = constantShift*ones(gridSize);
            obj = obj.cutAperture();
        end
    end
end