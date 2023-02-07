function [Summary, MitoObjectsGroupedPerWell] = f_imageAnalysisMito(ch1, ch2, ch3, WellThis, FieldThis, MesFile, PreviewPath, Layout)
%Summary is fine do previews and clean comments
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    Summary = table();

    %% segment nuclei
   
    NucLP = imfilter(ch1, fspecial('gaussian', 15, 3), 'symmetric'); % imtool(NucLP,[])
    %%figure; surf(fspecial('gaussian', 55, 11))
    NucMask = NucLP > 1000; % imtool(NucMask,[])
    NucMask = bwareaopen(NucMask, 750);
    [NucLM, NucCount] = bwlabeln(NucMask); % imtool(NucLM,[])
    NucObjects = regionprops('table', NucLM, ch1, {'MeanIntensity','Area','MajorAxisLength','MinorAxisLength','Perimeter'});
    NucSummary = table();
    NucSummary = varfun(@(x) mean(x), NucObjects);
    NucSummary.SumArea = sum(NucMask(:));
    NucSummary.Count = NucCount;
    NucSummary.Properties.VariableNames = strcat('Nuc_', NucSummary.Properties.VariableNames);

    
    %% Segment Lamp1
    Lamp1DoGKernel = fspecial('gaussian', 15, 1) - fspecial('gaussian', 15, 5);
    %figure; surf(Lamp1DoGKernel)% imtool(ch2,[])
    Lamp1DoG = imfilter(ch2, Lamp1DoGKernel,'symmetric');
    % imtool(Lamp1DoG, [])
    Lamp1Mask = Lamp1DoG>15 ;    
    Lamp1Mask = bwareaopen(Lamp1Mask, 5);
    %imtool(Lamp1Mask,[])
    [Lamp1LM, Lamp1Countall] = bwlabeln(Lamp1Mask);
    Lamp1ObjectsAll = regionprops('table', Lamp1LM, ch3, {'MeanIntensity','Area','MajorAxisLength','MinorAxisLength','Perimeter'});
    Lamp1SummaryAll = table();
    Lamp1SummaryAll = varfun(@(x) mean(x), Lamp1ObjectsAll);
    Lamp1SummaryAll.SumArea = sum(Lamp1Mask(:));
    Lamp1SummaryAll.Count = Lamp1Countall;
    Lamp1SummaryAll.Properties.VariableNames = strcat('Lamp1All_', Lamp1SummaryAll.Properties.VariableNames);
       

    %% Mitochondrial segmentation   
   
% ch1 = Dapi 
% ch2 = DeepRed AF647 Lamp1
% ch3 = Red AF568 Tom20

        
   DoGKernel = fspecial('gaussian', 15, 1) - fspecial('gaussian', 15, 5);
   %figure; surf(DoGKernel)
   MitoDoG = imfilter(ch3, DoGKernel,'symmetric') ;  %Tom20 is 568 
   % imtool(MitoDoG, []) ; 
   Mitomask = MitoDoG>150 ;    
   Mitomask = bwareaopen(Mitomask, 5);
   % imtool(Mitomask, [])
   
   Mitomatrix = bwlabeln(Mitomask);
   % imtool(Mitomatrix, []) 
    
   objects = regionprops('table', NucMask, ch3, {'Area','MeanIntensity'} );
   
    MitoMaskall= MitoDoG > 7; % imtool(MitoMaskall,[])
    MitoMaskall = MitoMaskall & ch3 > 0;
    MitoMaskall = bwareaopen(MitoMaskall, 5); 
   
    MitoMaskbright = MitoDoG > 10; % imtool(MitoMaskbright,[])
    MitoMaskbright = MitoMaskbright & ch3 > 0;
    MitoMaskbright = bwareaopen(MitoMaskbright, 5);
    [Mitoall, MitoCountall] = bwlabeln(MitoMaskall); %imtool(Mitoall,[])
    [Mitobright, MitoCountbright] = bwlabeln(MitoMaskbright); %imtool(Mitobright,[])
    MitoObjectsAll = regionprops('table', Mitoall, ch2, {'MeanIntensity','Area','MajorAxisLength','MinorAxisLength','Perimeter'});
    MitoObjectsAll.Properties.VariableNames = strcat('All_', MitoObjectsAll.Properties.VariableNames);
    MitoObjectsBright = regionprops('table', Mitobright, ch2, {'MeanIntensity','Area','MajorAxisLength','MinorAxisLength','Perimeter'});
    MitoObjectsBright.Properties.VariableNames = strcat('Bright_', MitoObjectsBright.Properties.VariableNames);
%     MitoSummary = table();
%     MitoSummary = varfun(@(x) mean(x), MitoObjects);
%     MitoSummary.SumArea = sum(MitoMask(:));
%     MitoSummary.Count = MitoCount;
%     MitoSummary.Properties.VariableNames = strcat('Mito_', MitoSummary.Properties.VariableNames);
    MitoSummaryAll = table();
    MitoSummaryAll = varfun(@(x) mean(x), MitoObjectsAll);
    MitoSummaryAll.SumArea = sum(MitoMaskall(:));
    MitoSummaryAll.Count = MitoCountall;
    MitoSummaryAll.Properties.VariableNames = strcat('MitoAll_', MitoSummaryAll.Properties.VariableNames);
    
    MitoSummaryBright = table();
    MitoSummaryBright = varfun(@(x) mean(x), MitoObjectsBright);
    MitoSummaryBright.SumArea = sum(MitoMaskbright(:));
    MitoSummaryBright.Count = MitoCountbright;
    MitoSummaryBright.Properties.VariableNames = strcat('MitoBright_', MitoSummaryBright.Properties.VariableNames);
   
  
   

         
   
    %% Morphmetric analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Skeletonize mitochondrial
    MitoMask = MitoMaskall; 
    Skeleton = Skeleton3D(MitoMask); %imtool(MitoMask + Skeleton, [])
    Branchpoints = bwmorph(Skeleton, 'branchpoints'); %imtool(MitoMask + Skeleton + Branchpoints, [])
    Endpoints = bwmorph(Skeleton, 'endpoints'); %imtool(MitoMask + Skeleton + Branchpoints + (2*Endpoints), [])


    %% Erode mitochondria for the bottleneck shape analysis

    MitoMaskEroded = imerode(MitoMask, strel('disk', 1));
    MitoMaskPerim = MitoMask - MitoMaskEroded; %imtool(MitoMaskPerim, []);

    %% Single mitochondria analysis

    [MitoStencil, MitoNumel] = bwlabel(MitoMask);% imtool(MitoStencil, [])
    MitoObjects = regionprops('table', MitoStencil, MitoStencil, {'MaxIntensity', 'BoundingBox', 'SubarrayIdx', 'PixelList', 'MajorAxisLength', 'MinorAxisLength'});
    MitoObjects.Properties.VariableNames = strrep(MitoObjects.Properties.VariableNames, 'MaxIntensity' ,'MitoIDinImage');
%     MitoObjects = regionprops('table', MitoStencil, CellsStencil, {'MaxIntensity', 'BoundingBox', 'SubarrayIdx', 'PixelList', 'MajorAxisLength', 'MinorAxisLength'});
%     MitoObjects.Properties.VariableNames{find(strcmp('MaxIntensity', MitoObjects.Properties.VariableNames))} = 'CellID';

    Area = {};
    Perimeter = {};
    FormFactor = {};
    AspectRatio = {};
    Skel = {};
    MitoObjectsNodeCount = {};
    MitoObjectsNodeDegree = {};
    MitoObjectsEndpointsCount = {};
    MitoShapeByPerimeter = {};
    MitoBodies = {};
    Tom20sumMito = {};

    for k = 1:height(MitoObjects)
        ConnectedPixels = [];
        ThisMitoRows = MitoObjects.SubarrayIdx{k,1};
        ThisMitoCols = MitoObjects.SubarrayIdx{k,2};
        ThisObjectMask = zeros([size(ThisMitoRows,2), size(ThisMitoCols,2)], 'logical');
        %Extract the relevant pixel indices within the bounding box
        PixelListInBB = MitoObjects.PixelList{k};
        PixelListRowsInBB = (PixelListInBB(:,2) - min(MitoObjects.SubarrayIdx{k,1}))+1;
        PixelListColsInBB = (PixelListInBB(:,1) - min(MitoObjects.SubarrayIdx{k,2}))+1;
        PixelListInBB = sub2ind(size(ThisObjectMask), PixelListRowsInBB, PixelListColsInBB);
        ThisObjectMask(PixelListInBB) = 1;
        %figure; imshow(imresize(ThisObjectMask,50, 'nearest'), [])

        ThisObjectSkel = Skeleton(ThisMitoRows, ThisMitoCols) .* ThisObjectMask;
        ThisObjectBranchpoints = Branchpoints(ThisMitoRows, ThisMitoCols) .* ThisObjectMask;
        ThisObjectEndpoints = Endpoints(ThisMitoRows, ThisMitoCols) .* ThisObjectMask;
        ThisObjectPerimeter = MitoMaskPerim(ThisMitoRows, ThisMitoCols) .* ThisObjectMask;
        ThisObjectBody = MitoMaskEroded(ThisMitoRows, ThisMitoCols) .* ThisObjectMask;
        ThisObjectTom20 = ch3(ThisMitoRows, ThisMitoCols) .* uint16(ThisObjectMask); % imtool(ThisObjectTom20, [])

        Area{k,1} = sum(ThisObjectMask(:));
        Perimeter{k,1} = sum(ThisObjectPerimeter(:));
        FormFactor{k,1} = (Perimeter{k,1}^2) / (4*pi*Area{k,1}); % AKA roundness or circularity: Lower values represent roundish objects
        AspectRatio{k,1} = MitoObjects.MajorAxisLength(k) / MitoObjects.MinorAxisLength(k);
        MitoShapeByPerimeter{k,1} = sum(ThisObjectBody(:)) / Perimeter{k,1};
        [~, MitoBodies{k,1}] = bwlabel(ThisObjectBody);
        Skel{k,1} = sum(ThisObjectSkel(:));
        Tom20sumMito{k,1} = sum(ThisObjectTom20(:));

        %Remove all Branchpoints and replace one by one to count directly
        %connected pixels
        ThisObjectSkelMinusBP = ThisObjectSkel - ThisObjectBranchpoints;
        %figure; imshow(imresize(ThisObjectSkelMinusBP,50, 'nearest'), [])
        BranchpointsCount = sum(ThisObjectBranchpoints(:));
        MitoObjectsNodeCount{k,1} = sum(ThisObjectBranchpoints(:));
        MitoObjectsEndpointsCount{k,1} = sum(ThisObjectEndpoints(:));

        if BranchpointsCount == 0
            MitoObjectsNodeDegree{k,1} = 0;
        else
            BP_Objects = regionprops('table', ThisObjectBranchpoints, {'PixelIdxList'});
            BP_IdxVec = BP_Objects.PixelIdxList;
            if iscell(BP_IdxVec); BP_IdxVec = BP_IdxVec{:}'; end
            for b = 1:size(BP_IdxVec,2)
                ThisObjectSkelMinusBPplusOne = ThisObjectSkelMinusBP;
                ThisObjectSkelMinusBPplusOne(BP_IdxVec(b)) = 3;
                %figure; imshow(imresize(ThisObjectSkelMinusBPplusOne,50, 'nearest'), [])
                ConnectionThisBP = regionprops('table', logical(ThisObjectSkelMinusBPplusOne), ThisObjectSkelMinusBPplusOne, {'PixelIdxList', 'MaxIntensity'});
                ConnectionThisBP = ConnectionThisBP(ConnectionThisBP.MaxIntensity == 3, :);
                ConnectedPixels(b) = size(ConnectionThisBP.PixelIdxList{:}, 1) - 1; % The BranchpointPixel is not counted
            end
            NodeDegree = mean(ConnectedPixels);
            MitoObjectsNodeDegree{k,1} = NodeDegree;
        end
    end
    %%
    MitoObjects = [MitoObjects, Tom20sumMito, Area, Perimeter, FormFactor, AspectRatio, MitoShapeByPerimeter, MitoBodies, Skel, MitoObjectsNodeCount, MitoObjectsEndpointsCount, MitoObjectsNodeDegree];
    if height(MitoObjects) == 0 % Skip well if no mitochondria were detected
        ExceptionIdx = ExceptionIdx + 1;
        Exceptions{ExceptionIdx} = [Barcode,'_',num2str(row),'_',num2str(column)];
        return
    end
    MitoObjects.Properties.VariableNames(7:end) = {'Tom20_MitoMask', 'MitoArea', 'MitoPerimeter', 'MitoFormFactor', 'MitoAspectRatio', 'MitoShapeByPerimeter', 'MitoErosionBodies', 'MitoSkel', 'MitoNodes', 'MitoEndpoints', 'MitoNodeDegree'};
    
    %% prepare for compact export of key features
    MitoObjects = MitoObjects(:, {'MitoIDinImage', 'MajorAxisLength', 'MinorAxisLength', 'Tom20_MitoMask', 'MitoArea', 'MitoPerimeter', 'MitoFormFactor', 'MitoAspectRatio',...
                                  'MitoShapeByPerimeter', 'MitoErosionBodies', 'MitoSkel', 'MitoNodes', 'MitoEndpoints', 'MitoNodeDegree'});
                              
    Well2table = cell(height(MitoObjects), 1); Well2table = cellfun(@(x) WellThis, Well2table, 'UniformOutput', false); 
    Field2table = cell(height(MitoObjects), 1); Field2table = cellfun(@(x) FieldThis, Field2table, 'UniformOutput', false);
    %Barcode2table = cell(height(MitoObjects), 1); Barcode2table = cellfun(@(x) MetaThis.Barcode, Barcode2table, 'UniformOutput', false);
    %CellLine2table = cell(height(MitoObjects), 1); CellLine2table = cellfun(@(x) MetaThis.CellLine, CellLine2table, 'UniformOutput', false);
    %Condition2table = cell(height(MitoObjects), 1); Condition2table = cellfun(@(x) MetaThis.Barcode, Condition2table, 'UniformOutput', false);
    
    
    MitoObjects.Well = Well2table;
    MitoObjects.Field = Field2table;
    %MitoObjects.Barcode = Barcode2table;
    %MitoObjects.CellLine = CellLine2table;
    %MitoObjects.Condition = Condition2table;


    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Collect features    
    
    MetaColumns = table();
    MetaColumns.Well = WellThis;
    MetaColumns.Field = FieldThis;
    MetaColumns.Mes = {MesFile};
%     MetaColumns.Barcode = {MetaThis.Barcode};
%     MetaColumns.CellLine = {MetaThis.CellLine};
%     MetaColumns.Condition = {MetaThis.ExperimentalCondition};
%     
    

    %% Collect outputs
    %Summary = [MetaColumns, MitoSummaryAll, MitoSummaryBright, NucSummary];
    Summary = [MetaColumns, MitoSummaryAll, Lamp1SummaryAll, MitoSummaryBright, NucSummary];
    Summary.Well = {Summary.Well};
    MitoObjectsGroupedPerWell = grpstats(MitoObjects(:, [2:15]), {'Well'}, {'mean', 'std'});
    Summary = outerjoin(Summary, MitoObjectsGroupedPerWell);

    %% Previews
    
    % Name
    LayoutThis = Layout(strcmp(Layout.Well, WellThis),:);
    NameThis =  [LayoutThis.CellLine{:}, '__', LayoutThis.ExperimentalCondition{:}];

% ch2 = DeepRed AF647 Lamp1
% ch3 = Red AF568 Tom20

    imSize = size(MitoMaskall);
    [BarMask, BarCenter] = f_barMask(20, 0.32393102760889064, imSize, imSize(1)-50, 50, 20);
    %imtool(BarMask)

    NucPreview = f_imoverlayIris(imadjust(ch1, [0 0.01], [0 1]), imdilate(bwperim(NucMask),strel('disk', 1)), [0 0 1]);
    NucPreview = f_imoverlayIris(NucPreview, BarMask, [1 1 1]);
    %imtool(NucPreview)
    
    Lamp1Preview = f_imoverlayIris(imadjust(ch2, [0 0.015], [0 1]), bwperim(Lamp1Mask), [0 1 0]);
    Lamp1Preview = f_imoverlayIris(Lamp1Preview, BarMask, [1 1 1]);
    %imtool(Lamp1Preview)


    MitoPreviewTom20All = f_imoverlayIris(imadjust(ch3, [0 0.03], [0 1]), bwperim(MitoMaskall), [1 0 0]);
    MitoPreviewTom20All = f_imoverlayIris(MitoPreviewTom20All, BarMask, [1 1 1]);
    %imtool(MitoPreviewTom20All)

    %MitoPreviewAll = f_imoverlayIris(imadjust(ch2, [0 0.03], [0 1]), bwperim(MitoMask), [1 0 0]);
    %MitoPreviewAll = f_imoverlayIris(imadjust(ch2, [0 0.005], [0 1]), bwperim(MitoMaskall), [1 0 0]);
    MitoPreviewAll = f_imoverlayIris(imadjust(ch2, [0 0.15], [0 1]), bwperim(MitoMaskall), [1 0 0]);
    MitoPreviewAll = f_imoverlayIris(MitoPreviewAll, BarMask, [1 1 1]);
    %imtool(MitoPreviewAll)
    %imtool(ch2, []) % Lamp1
    

    MitoPreviewTom20Bright = f_imoverlayIris(imadjust(ch3, [0 0.03], [0 1]), bwperim(MitoMaskbright), [1 0 0]);
    MitoPreviewTom20Bright = f_imoverlayIris(MitoPreviewTom20Bright, BarMask, [1 1 1]);
    %imtool(MitoPreviewTom20Bright)

    %MitoPreviewBright = f_imoverlayIris(imadjust(ch2, [0 0.03], [0 1]), bwperim(MitoMask), [1 0 0]);
    MitoPreviewBright = f_imoverlayIris(imadjust(ch2, [0 0.05], [0 1]), bwperim(MitoMaskbright), [1 0 0]);
    MitoPreviewBright = f_imoverlayIris(MitoPreviewBright, BarMask, [1 1 1]);
    %imtool(MitoPreviewBright)
    
    % Skeleton
    PreviewSkeleton = imoverlay2(MitoMask, Skeleton, [0 0 1]);
    PreviewSkeleton = imoverlay2(PreviewSkeleton, Branchpoints, [0 1 0]);
    PreviewSkeleton = imoverlay2(PreviewSkeleton, Endpoints, [1 0 0]);
    PreviewSkeleton = imoverlay2(PreviewSkeleton, BarMask, [1 1 1]); %imtool(PreviewSkeleton, [])
    
    
    
%     RGBPreviewMito = cat(3, imadjust(medfilt3(ch3),[0.002 0.065], [0 1]), imadjust(medfilt3(ch2),[0.01 0.02], [0 1]),imadjust(medfilt3(ch1),[0.001 0.15], [0 1]) );
%     RGBPreviewMito = imoverlay(RGBPreviewMito, BarMask, [1 1 1]); % imtool(RGBPreviewMito)
   

    RGBPreviewMito = cat(3, imadjust(medfilt3(ch3),[0.005 0.015], [0 1]), imadjust(medfilt3(ch2),[0.005 0.013], [0 1]),imadjust(medfilt3(ch1),[0.01 0.04], [0 1]) );
    RGBPreviewMito = imoverlay(RGBPreviewMito, BarMask, [1 1 1]); % imtool(RGBPreviewMito)
    
    RGBPreviewMitoImad = cat(3, imadjust(medfilt3(ch3)), imadjust(medfilt3(ch2)), imadjust(medfilt3(ch1)) );
    RGBPreviewMitoImad = imoverlay(RGBPreviewMitoImad, BarMask, [1 1 1]); % imtool(RGBPreviewMitoImad)
    
    
    RGBPreviewMitoPathImad = [PreviewPath, filesep, WellThis, '_', FieldThis, '_', NameThis, '_RGBMito_imadjust.png'];
    RGBPreviewMitoPath = [PreviewPath, filesep, WellThis, '_', FieldThis, '_', NameThis, '_RGBMito.png'];
    NucPreviewPath = [PreviewPath, filesep, WellThis, '_', FieldThis, '_', NameThis, '_Nuc.png'];
    MitoLamp1PreviewPathAll = [PreviewPath, filesep, WellThis, '_', FieldThis, '_', NameThis, '_MitoTom20All.png'];
    MitoLamp1PreviewPathBright = [PreviewPath, filesep, WellThis, '_', FieldThis, '_', NameThis, '_MitoTom20Bright.png'];
    MitoPreviewPathAll = [PreviewPath, filesep, WellThis, '_', FieldThis, '_', NameThis, '_MitoAll.png'];
    MitoPreviewPathBright = [PreviewPath, filesep, WellThis, '_', FieldThis, '_', NameThis, '_MitoBright.png'];
    MitoSkeletonPreviewPath = [PreviewPath, filesep, WellThis, '_', FieldThis, '_', NameThis, '_Skeleton.png'];
    Lamp1PreviewPath = [PreviewPath, filesep, WellThis, '_', FieldThis, '_', NameThis, '_Lamp1.png'];

    
    imwrite(RGBPreviewMitoImad, RGBPreviewMitoPathImad)
    imwrite(RGBPreviewMito, RGBPreviewMitoPath)
    imwrite(NucPreview, NucPreviewPath)
    imwrite(MitoPreviewTom20All, MitoLamp1PreviewPathAll)
    imwrite(MitoPreviewTom20Bright, MitoLamp1PreviewPathBright)
    imwrite(MitoPreviewAll, MitoPreviewPathAll)
    imwrite(MitoPreviewBright, MitoPreviewPathBright)
    imwrite(PreviewSkeleton, MitoSkeletonPreviewPath);
    imwrite(Lamp1Preview, Lamp1PreviewPath);


end


     

