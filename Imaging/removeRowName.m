function [x] = removeRowName(x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if istable(x)
       x.Properties.RowNames = {};    
    end
end

