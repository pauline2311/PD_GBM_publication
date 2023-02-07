function tableAnnotated = Iris_AnnotateTable(tableObjects, LayoutFromDB, ColumnsOrAll)
% Append metadata columns to a table with well labels in the style A01...H12
% author: Paul Antony 20210112
%   ColumnsOrAll: set 'All' to append all columns
%                 set cellArray of features to append specific columns
%                 existing in LayoutFromDB

if strcmp(ColumnsOrAll, 'All')
    Columns = LayoutFromDB.Properties.VariableNames;
else
    Columns = ColumnsOrAll;
end

if ~ismember('Well', Columns)
    Columns = [Columns, {'Well'}];
end

tableAnnotated = outerjoin(tableObjects, LayoutFromDB(:, Columns), 'Type', 'left', 'MergeKeys', true);

