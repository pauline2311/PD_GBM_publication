function [Objects] = Iris_ObjectsInsertSubObjects(StencilMainObjects, MaskSubObjects, connectivity, varargin)
%Insert images of submasks and basic sub-mask morphology features into an Objectlist containing {'PixelIdxList'}
% Paul Antony 20210203
%   connectivity (4 or 8)


    Objects = regionprops('table', StencilMainObjects, {'Image', 'PixelIdxList', 'PixelList', 'SubarrayIdx'});
    if height(Objects) == 0
        Objects = table();
        return
    end
    
    Objects.SubMasks = rowfun(@(x,y) MaskSubObjects(x{:}) .* y{:}, Objects, 'InputVariables', {'SubarrayIdx', 'Image'}, 'OutputFormat', 'cell');  
    
    %% Insert basic morphology features
    Objects.SubArea = rowfun(@(x) sum(x{:}(:)), Objects, 'InputVariables', {'SubMasks'}, 'OutputFormat', 'cell');
    Objects.CountSubObjects = rowfun(@(x) Iris_CountObjects(x{:}, 8), Objects, 'InputVariables', {'SubMasks'}, 'OutputFormat', 'cell');

    if nargin -3 > 0
        for i = 1:nargin-3
            disp('Intensity analysis for optional channel ...')
            OptChThis = varargin{i};
            OptionalChThisIntensityInSubMask = rowfun(@(x,y) uint16(OptChThis(x{:})) .* uint16(y{:}), Objects, 'InputVariables', {'SubarrayIdx', 'SubMasks'}, 'OutputFormat', 'cell');  
            OptionalChThisIntensityMeanInSubMask = cellfun(@(a) mean(a(a>0)), OptionalChThisIntensityInSubMask);
            eval(['Objects.Varargin', num2str(i), '_meanSubIntensity = OptionalChThisIntensityMeanInSubMask;']);
        end
    end

end

