function [Layout] = Iris_GetLayout(InPath)
%Load layout saved in hcsdb which has been uploaded via IrisStart.py
%   author: paul.antony@uni.lu 20210107
%   InPath: the path corresponding to InPathIris in hcsdb
    Layout = readtable([InPath, filesep, 'LayoutHCSDB.csv'], 'Delimiter', ',');
    try
        Layout = Layout(cellfun(@(x) ~isempty(x), Layout.ExperimentalCondition),:);
    catch % no experimental condition used, example only cell-lines
        disp('No Experimental condition defined')
    end
    Layout.Well = rowfun(@(r,c) sprintf('%s%02d', r, c), Layout, 'Inputvariables', {'RowLetter', 'ColumnNumber'}, 'ExtractCellContents', true, 'OutputFormat', 'cell')
end

