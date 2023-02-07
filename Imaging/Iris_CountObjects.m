function [Count] = Iris_CountObjects(Mask,connectivity)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    [~, Count] = bwlabeln(Mask, connectivity);
end

