function [CellObjects] = f_imageAnalysisAstro(ch1, ch2, ch3, WellThis, FieldThis, MesFile, PreviewPath, Layout)
%Summary is fine do previews and clean comments
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    Summary = table();

    %% segment nuclei
    
    %ch1 = max(ch1all, [], 3); % imtool(ch1, [])
    %ch1 = ch1all(:,:,t); 
    NucLP = imfilter(ch1, fspecial('gaussian', 15, 3), 'symmetric'); % imtool(NucLP,[])
    %%figure; surf(fspecial('gaussian', 55, 11))
    NucMask = NucLP > 1000; % imtool(NucMask,[])
    NucMask = bwareaopen(NucMask, 750);
    
     %% Split Nuclei
     D = bwdist(~NucMask);
     %imtool(D,[]);
     WI = imcomplement(D);
     W = watershed(imhmin(WI, 1));%it(W)
     NucStencil = uint16(NucMask) .* uint16(W);
     NucMask = NucStencil > 0;
     NucMask = bwareaopen(NucMask, 750);


     [NucLM, NucCount] = bwlabeln(NucMask, 4); %imtool(NucLM,[])
     NucObjects = regionprops('table', NucLM, ch1, {'MeanIntensity','Area','MajorAxisLength','MinorAxisLength','Perimeter'});
     NucSummary = table();
     NucSummary = varfun(@(x) mean(x), NucObjects);
     NucSummary.SumArea = sum(NucMask(:));
     NucSummary.Count = NucCount;
     NucSummary.Properties.VariableNames = strcat('Nuc_', NucSummary.Properties.VariableNames);

    
    %% Perinuc analysis    
    %% Dilate and resplit nuclei
    D = bwdist(NucMask);% imtool(D, [])
    W = watershed(D);% imtool(W, [])
    NucStencil = uint16(NucMask) .* uint16(W); % imtool(NucStencil, [])
    NucExtended = uint16(imdilate(NucMask, strel('disk', 15))) .* uint16(W); % imtool(NucExtended, [])
    Perinuc = NucExtended .* uint16(~NucMask);% imtool(Perinuc, [])
    PerinucMask = logical(Perinuc);
    PerinucMask = bwareaopen(PerinucMask, 500);% imtool(PerinucMask)
    PerinucLM = bwlabeln(PerinucMask);
    
    
    %% Intensities by channel 
    % CellMask
    CellMask = NucMask | PerinucMask; %imtool(CellMask,[])
    CellStencil = bwlabeln(CellMask); %imtool(CellStencil,[])
    CellObjectsNuc = Iris_ObjectsInsertSubObjects(CellStencil, NucMask, 4);
    CellObjectsNuc.Properties.VariableNames = strrep(CellObjectsNuc.Properties.VariableNames, 'SubMasks','SubMasksNuc');
    CellObjectsPeriNuc = Iris_ObjectsInsertSubObjects(CellStencil, PerinucMask, 4);
    CellObjectsPeriNuc.Properties.VariableNames = strrep(CellObjectsPeriNuc.Properties.VariableNames, 'SubMasks','SubMasksPeriNuc');
    CellObjects = CellObjectsNuc;
    CellObjects.SubMasksPeriNuc = CellObjectsPeriNuc.SubMasksPeriNuc;
    %CellObjects.SubMaskCell = rowfun(@(A, B) {A{:}|B{:}} , CellObjects, 'InputVariables', {'SubMasksNuc', 'SubMasksPeriNuc'}, 'OutputFormat', 'cell', 'ExtractCellContents', true);
    CellObjects.SubMaskCell = rowfun(@(A, B) A|B , CellObjects, 'InputVariables', {'SubMasksNuc', 'SubMasksPeriNuc'}, 'OutputFormat', 'cell', 'ExtractCellContents', true);
    
    NucRedThis = [];
    PerinucRedThis = [];
    CellRedThis = [];
    NucDeepRedThis = [];
    PerinucDeepRedThis = [];
    CellDeepRedThis = [];
    
    PeriNucNucDeepRedThis = [];
    PeriNucNucRedThis = [];
    
    
    % imtool(CellObjects{1, 'Image'}{:}, [])
    for i = 1:height(CellObjects)
        BlueThis = ch1(CellObjects{i,'SubarrayIdx'}{1}, CellObjects{i,'SubarrayIdx'}{2}); %imtool(BlueThis, []) 
        RedThis = ch2(CellObjects{i,'SubarrayIdx'}{1}, CellObjects{i,'SubarrayIdx'}{2}); %imtool(RedThis, [])
        DeepRedThis = ch3(CellObjects{i,'SubarrayIdx'}{1}, CellObjects{i,'SubarrayIdx'}{2}); %imtool(DeepRedThis, [])
        NucRedThis(i) = mean(RedThis(find(CellObjects{i,'SubMasksNuc'}{:}))); % imtool(NucRedThis,[])
        PerinucRedThis(i) = mean(RedThis(find(CellObjects{i,'SubMasksPeriNuc'}{:}))); % imtool(NucRedThis,[])
        CellRedThis(i) = mean(RedThis(find(CellObjects{i,'Image'}{:}))); % imtool(NucRedThis,[])
        NucDeepRedThis(i) = mean(DeepRedThis(find(CellObjects{i,'SubMasksNuc'}{:}))); % imtool(NucRedThis,[])
        PerinucDeepRedThis(i) = mean(DeepRedThis(find(CellObjects{i,'SubMasksPeriNuc'}{:}))); % imtool(NucRedThis,[])
        CellDeepRedThis(i) = mean(DeepRedThis(find(CellObjects{i,'Image'}{:}))); % imtool(NucRedThis,[])
        
        
        PeriNucNucDeepRedThis(i) = mean(DeepRedThis(find(CellObjects{i,'SubMaskCell'}{:}))); % imtool(PeriNucNucDeepRedThis,[])
        PeriNucNucRedThis(i) = mean(RedThis(find(CellObjects{i,'SubMaskCell'}{:}))); % imtool(PeriNucNucRedThis,[])
        
               
        %NucRedThisbyNuc(i) = sum(NucRedThis(:)) / sum(BlueThis(:)); % Sum of pixels in Mask normalized to Nuc
        
        
        
        
    end
    
    CellObjects = CellObjects(:,'CountSubObjects');
    CellObjects.Well = repmat({WellThis},[height(CellObjects), 1]);
    CellObjects = Iris_AnnotateTable(CellObjects, Layout, 'All');
    CellObjects.NucArea = repmat(sum(NucMask(:)),[height(CellObjects), 1]);
    CellObjects.NucRed = NucRedThis';
    CellObjects.PerinucRed = PerinucRedThis';
    CellObjects.CellRed = CellRedThis';
    CellObjects.NucDeepRed = NucDeepRedThis';
    CellObjects.PerinucDeepRed = PerinucDeepRedThis';
    CellObjects.CellDeepRed = CellDeepRedThis';
    
    CellObjects.PeriNucNucDeepRed = PeriNucNucDeepRedThis';
    CellObjects.PeriNucNucRed = PeriNucNucRedThis';
    
    %CellObjects.NucRedbyNuc = NucRedThisbyNuc';
    
    
    
    %% Previews    
    % Name
    LayoutThis = Layout(strcmp(Layout.Well, WellThis),:);
    NameThis =  [LayoutThis.CellLine{:}, '__', LayoutThis.ExperimentalCondition{:}];
    
   
    
    imSize = size(ch1);
    [BarMask, BarCenter] = f_barMask(15, 0.10758027143330025, imSize, imSize(1)-50, 50, 10);
    %it(BarMask)
    
      
%     RGBPreview = cat(3, imadjust(ch3, [0 0.02], [0 1]), imadjust(ch2, [0 0.02], [0 1]), imadjust(ch1, [0 0.01], [0 1]));
%     RGBPreview = cat(3, imadjust(ch3, [0 0.1], [0 1]), imadjust(ch2, [0 0.1], [0 1]), zeros(size(ch1), 'uint16'));
%     RGBPreview = f_imoverlayIris(RGBPreview, BarMask, [1 1 1]);
%     %imtool(RGBPreview)
%      
    
    RGBPreview = cat(3, imadjust(medfilt3(ch3),[0.002 0.065], [0 1]), imadjust(medfilt3(ch2),[0.01 0.02], [0 1]),imadjust(medfilt3(ch1),[0.001 0.15], [0 1]) );
    RGBPreview = imoverlay(RGBPreview, BarMask, [1 1 1]); % imtool(RGBPreview)
    
    RGBPreviewPM = cat(3, imadjust(medfilt3(ch3),[0.005 0.015], [0 1]), imadjust(medfilt3(ch2),[0.005 0.013], [0 1]),imadjust(medfilt3(ch1),[0.01 0.04], [0 1]) );
    RGBPreviewPM = imoverlay(RGBPreviewPM, BarMask, [1 1 1]); % imtool(RGBPreviewPM)   
    
    RGBPreviewImad = cat(3, imadjust(ch3), imadjust(ch2), imadjust(ch1) );
    RGBPreviewImad = imoverlay(RGBPreviewImad, BarMask, [1 1 1]); % imtool(RGBPreviewImad)
   
    
   
    NucPreview = f_imoverlayIris(imadjust(ch1, [0 0.1], [0 1]), imdilate(bwperim(NucMask),strel('disk', 1)), [0 0 1]);
    NucPreview = f_imoverlayIris(NucPreview, BarMask, [1 1 1]);
    %imtool(NucPreview)
    
    
    
    RedPreviewAuto = f_imoverlayIris(imadjust(ch3), imdilate(bwperim(PerinucMask),strel('disk', 1)), [0 0 1]);
    RedPreviewAuto = f_imoverlayIris(RedPreviewAuto, BarMask, [1 1 1]); %imtool(RedPreviewAuto)
    
    
    RedPreview = f_imoverlayIris(imadjust(ch3, [0 0.1], [0 1]), imdilate(bwperim(PerinucMask),strel('disk', 1)), [0 0 1]);
    %RedPreview = f_imoverlayIris(RedPreview, PerinucMask, [0 1 1]);
    RedPreview = f_imoverlayIris(RedPreview, BarMask, [1 1 1]); %imtool(RedPreview)
    
    
    DeepRedPreviewAuto = f_imoverlayIris(imadjust(ch2), imdilate(bwperim(PerinucMask),strel('disk', 1)), [0 0 1]);
    DeepRedPreviewAuto = f_imoverlayIris(DeepRedPreviewAuto, BarMask, [1 1 1]); %imtool(DeepRedPreviewAuto)
    
    
    DeepRedPreview = f_imoverlayIris(imadjust(ch2, [0 0.1], [0 1]), imdilate(bwperim(PerinucMask),strel('disk', 1)), [0 0 1]);
    %DeepRedPreview = f_imoverlayIris(DeepRedPreview, PerinucMask, [0 1 1]);
    DeepRedPreview = f_imoverlayIris(DeepRedPreview, BarMask, [1 1 1]); %imtool(DeepRedPreview)
    
    
    
    %NucPreviewPath = [PreviewPath, filesep, '_', FieldThis, '_Nuc.png'];
    %RGBPreviewPath = [PreviewPath, filesep, '_', FieldThis, '_RGB.png'];
    
    %LayoutPlusThisWell = Layout{strcmp(Layout.Well, WellThis), 'ExperimentalCondition'}{:};

%     NucPreviewPath = [PreviewPath, filesep, LayoutPlusThisWell, '_', WellThis, '_', FieldThis, '_', NameThis, '_Nuc.png'];
%     RGBPreviewPath = [PreviewPath, filesep, LayoutPlusThisWell, '_', WellThis, '_', FieldThis, '_', NameThis, '_RGB.png'];

    NucPreviewPath = [PreviewPath, filesep, WellThis, '_', FieldThis, '_', NameThis, '_Nuc.png'];
    RGBPreviewPath = [PreviewPath, filesep, WellThis, '_', FieldThis, '_', NameThis, '_RGB.png'];
    RGBPreviewPMPath = [PreviewPath, filesep, WellThis, '_', FieldThis, '_', NameThis, '_RGB_PM.png'];
    RGBPreviewPathImad = [PreviewPath, filesep, WellThis, '_', FieldThis, '_', NameThis, '_RGB_imadjust.png'];
    RedPreviewPath = [PreviewPath, filesep, WellThis, '_', FieldThis, '_', NameThis, '_Red.png'];
    DeepRedPreviewPath = [PreviewPath, filesep, WellThis, '_', FieldThis, '_', NameThis, '_DeepRed.png'];
    RedPreviewAutoPath = [PreviewPath, filesep, WellThis, '_', FieldThis, '_', NameThis, '_RedAuto.png'];    
    DeepRedPreviewAutoPath = [PreviewPath, filesep, WellThis, '_', FieldThis, '_', NameThis, '_DeepRedAuto.png'];
    
    
    
    
    imwrite(NucPreview, NucPreviewPath)
    imwrite(RGBPreview, RGBPreviewPath)
    imwrite(RGBPreviewPM, RGBPreviewPMPath)
    imwrite(RGBPreviewImad, RGBPreviewPathImad)
    imwrite(RedPreview, RedPreviewPath)
    imwrite(DeepRedPreview, DeepRedPreviewPath)
    imwrite(RedPreviewAuto, RedPreviewAutoPath)
    imwrite(DeepRedPreviewAuto, DeepRedPreviewAutoPath)


end

