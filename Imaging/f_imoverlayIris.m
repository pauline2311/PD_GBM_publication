function OverlayRGB = f_imoverlayIris(Im, Mask, RGB)
%IMOVERLAY Adds an overlay mask to a color image without using interactive/graphical Matlab funcctionality.
%
%   Recommended for Slurm sbatch runs
%   Im has to be two dimensional
%   Mask has to be twoo dimensional as well
%   RGB is a red green blue vector with a range from 0 to 1 for example [1 0 0]
%   Example: Im = im2uint16(rand(100)); Mask = Im > 30000; OverlayRGB = f_imoverlayIris(Im, Mask, [1 0.5 0]); imtool(OverlayRGB)

    Mask = logical(Mask);
    ImClass = class(Im);
    switch ImClass
        case 'uint8'
            MaxPossible = (2^8)-1;
        case 'uint16'
            MaxPossible = (2^16)-1;
        otherwise
            disp('Please input the graytone image as uint8 or as uint16')
            return
    end
    if size(Im, 3) == 1
        ImR = Im;
        ImG = Im;
        ImB = Im;
    else
        ImR = Im(:,:,1);
        ImG = Im(:,:,2);
        ImB = Im(:,:,3);
    end
    ImR(Mask) = RGB(1) * MaxPossible;
    ImG(Mask) = RGB(2) * MaxPossible;
    ImB(Mask) = RGB(3) * MaxPossible;

    OverlayRGB = cat(3, ImR, ImG, ImB);



end