%% Collect Linux\Slurm metadata
disp('Paul can edit here')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Job is running on node:')
[~, node] = system('hostname');
disp(node)
disp('Job is run by user:')
[~, user] = system('whoami');
disp(user)
disp('Current slurm jobs of current user:')
[~, sq] = system(['squeue -u ', user]);
disp(sq)
tic
disp(['Start: ' datestr(now, 'yyyymmdd_HHMMSS')])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


%% Loading library 

addpath(genpath('/work/projects/lcsb_hcs/Library/hcsforge'))
addpath(genpath('/work/projects/lcsb_hcs/Library/hcsIris'))


%% Load Data 

if ~exist('InPath') % if Inpath is provided  via command line, use that one
    %InPath = '/work/projects/lcsb_hcs/Data/PaulineMencke/AstroMitoMorph/PlateTest_20210126_140544_in';
    InPath = '/work/projects/lcsb_hcs/Data/PaulineMencke/AstroMitoMorph/20210504_PM_AstroICC_IBO_N1_20210504_094756_in';
    
end

Layout = Iris_GetLayout(InPath);

MesPath = ls([InPath, '/*.mes']); MesPath = MesPath(1:end-1); % remove line break
MetaData = f_CV8000_getChannelInfo(InPath, MesPath);

if ~exist('OutPath') % if Outpath is provided  via command line, use that one
    %OutPath = '/work/projects/lcsb_hcs/Data/PaulineMencke/AstroMitoMorph/PlateTest_20210126_140544_out';
    OutPath = '/work/projects/lcsb_hcs/Data/PaulineMencke/AstroMitoMorph/20210504_PM_AstroICC_IBO_N1_20210504_094756_out2';
    
end

%% Prepare folders
mkdir(OutPath)
PreviewPath = [OutPath, filesep, 'Previews'];
mkdir(PreviewPath)

%% Log
f_LogDependenciesLinux(mfilename, OutPath)


%% Load Metadata
ObjectsAll = {};
MitoObjectsAll = {};
SummaryAll = {};
%MetaData = f_CV8000_getChannelInfo(InPath, MesPath);
InfoTable = MetaData.InfoTable{:};
Wells = unique(InfoTable.Well);
fieldProgress = 0;
for w = 1:numel(Wells)
    WellThis = Wells{w};
    InfoTableThisWell = InfoTable(strcmp(InfoTable.Well, WellThis),:);
    FieldsThisWell = unique(InfoTableThisWell.Field);
    for f = 1:numel(FieldsThisWell)
        fieldProgress = fieldProgress + 1;
        FieldThis = FieldsThisWell{f};
        InfoTableThisField = InfoTableThisWell(strcmp(InfoTableThisWell.Field, FieldsThisWell{f}),:);
        ChannelsThisField =  unique(InfoTableThisField.Channel);
        ImPaths = cell(1, numel(ChannelsThisField));
        for c = 1:numel(ChannelsThisField)
            ChannelThis = ChannelsThisField{c};
            InfoTableThisChannel = InfoTableThisField(strcmp(InfoTableThisField.Channel,ChannelThis),:);
            InfoTableThisChannel = sortrows(InfoTableThisChannel, 'Plane', 'ascend');
            chThisPaths = cell(numel(ChannelsThisField),1);
            for p = 1:height(InfoTableThisChannel)
                chThisPaths{p} = InfoTableThisChannel{p, 'file'}{:};
                %for t = 1:height()
            end
            ImPaths{c} = chThisPaths;
            MesFile = MetaData.MeasurementSettingFileName;
        end
       FieldMetaData{fieldProgress} = {ImPaths, MesFile, Wells{w}, FieldsThisWell{f}};
    end
end

disp('Debug point')
FieldMetaDataTable = cell2table(vertcat(FieldMetaData{:}))
LayoutPlus = Layout;
LayoutPlus.IsMitoAssay = rowfun(@(x) ~isempty(regexp(x, '.*Tom20.*', 'tokens')), LayoutPlus, 'InputVariables', 'ExperimentalCondition', 'ExtractCellContents', true, 'OutputFormat', 'cell')
%~isempty(regexp('Lamp1_647_Tom20_568_120minFCCP', '.*Tom20.*', 'tokens'))
LayoutPlus.IsNeuronsAssay = rowfun(@(x) ~isempty(regexp(x, '.*Neurons.*', 'tokens')), LayoutPlus, 'InputVariables', 'CellLine', 'ExtractCellContents', true, 'OutputFormat', 'cell')

%/work/projects/lcsb_hcs/Advanced/PaulineMencke/AstrocytePaper/ICC_AstroMitoMorph/
parfor i = 1:numel(FieldMetaData)

%for i = 3172
%for i = [1703, 2129]
%for i = 126
%for i = 2441
%for i = 229
%for i = 1489
%for i = 439
%for i=[3339, 3445:3455]
%for i = 1:numel(FieldMetaData)

    try
        i    
        ch1files = sort(FieldMetaData{i}{1}{1}(1));
        ch1Collector = cellfun(@(x) imread(x), ch1files, 'UniformOutput', false);
        ch1 = cat(3,ch1Collector{:}); % vol(ch1, 0, 2000) Hoechst
    
        ch2files = sort(FieldMetaData{i}{1}{2}(1));
        ch2Collector = cellfun(@(x) imread(x), ch2files, 'UniformOutput', false);
        ch2 = cat(3,ch2Collector{:}); % vol(ch2, 0, 800) DeepRed
        
        ch3files = sort(FieldMetaData{i}{1}{3}(1));
        ch3Collector = cellfun(@(x) imread(x), ch3files, 'UniformOutput', false);
        ch3 = cat(3,ch3Collector{:}); % vol(ch3, 0, 2000) Red

        MesFile = FieldMetaData{i}{2};
        WellThis = FieldMetaData{i}{3};
        FieldThis = FieldMetaData{i}{4};
        LayoutPlusThisWell = LayoutPlus(strcmp(LayoutPlus.Well, WellThis), :);
        
        
        if LayoutPlusThisWell.IsMitoAssay{:}
            %continue
            [Summary, MitoObjectsGroupedPerWell] = f_imageAnalysisMito(ch1, ch2, ch3, WellThis, FieldThis, MesFile, PreviewPath, Layout);
            MitoObjectsAll{i} = MitoObjectsGroupedPerWell;
            SummaryAll{i} = Summary;
%         elseif LayoutPlusThisWell.IsNeuronsAssay{:}            
%             ObjectsNeurons = f_imageAnalysisNeurons(ch1, ch2, ch3, WellThis, FieldThis, MesFile, PreviewPath, Layout);
%             ObjectsAllNeurons{i} = ObjectsNeurons;            
        else
            Objects = f_imageAnalysisAstro(ch1, ch2, ch3, WellThis, FieldThis, MesFile, PreviewPath, Layout);
            ObjectsAll{i} = Objects;
            
        end
        
        
        %ObjectsAll{i} = Objects;
    catch mygreaterror        
       disp(mygreaterror)       
       continue
    end
end

%%
save([OutPath, filesep, 'data.mat'], 'MitoObjectsAll', 'SummaryAll', 'ObjectsAll');


MitoObjectsAll = cellfun(@(x) removeRowName(x), MitoObjectsAll, 'UniformOutput', false);
MitoObjectsAll = vertcat(MitoObjectsAll{:});
writetable(MitoObjectsAll, [OutPath, filesep, 'MitoObjectsAll.csv'])

SummaryAll = cellfun(@(x) removeRowName(x), SummaryAll, 'UniformOutput', false);
SummaryAll = vertcat(SummaryAll{:});
writetable(SummaryAll, [OutPath, filesep, 'SummaryAll.csv'])

% ObjectsAllNeurons = cellfun(@(x) removeRowName(x), ObjectsAllNeurons, 'UniformOutput', false);
% ObjectsAllNeurons = vertcat(ObjectsAllNeurons{:});
% writetable(ObjectsAllNeurons, [OutPath, filesep, 'ObjectsAllNeurons.csv'])

ObjectsAll = cellfun(@(x) removeRowName(x), ObjectsAll, 'UniformOutput', false);
ObjectsAll = vertcat(ObjectsAll{:});
writetable(ObjectsAll, [OutPath, filesep, 'ObjectsAll.csv'])
disp('Script completed successfully')
