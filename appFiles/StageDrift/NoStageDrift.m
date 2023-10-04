classdef NoStageDrift < StageDrift

    methods
        function obj = NoStageDrift()
            obj.motion = Length([0 0 0],'nm');
        end
    end
end