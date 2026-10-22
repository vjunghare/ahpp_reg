% AHPP regression: predict the IC50 (uM) of peptides with a regression decision tree.
%
% The tree is trained every time the script runs, from the training data in
% data/train/data_regr2.mat, so no saved model is needed. It is then applied to the
% peptides in the two test files in data/example_test/.
%
% Requires MATLAB with the Statistics and Machine Learning Toolbox (fitrtree, predict).
% Run it from any folder; the file locations below are relative to this script.

clear; clc;
here = fileparts(mfilename('fullpath'));

%% Settings
trainFile   = fullfile(here, 'data', 'train', 'data_regr2.mat');
structFile  = fullfile(here, 'data', 'example_test', 'ml_csv.csv');             % 8 structure features per peptide
pseaacFile  = fullfile(here, 'data', 'example_test', 'peptides.pfpseaac.csv');  % 23 sequence features per peptide
outFile     = fullfile(here, 'ic50_peptide.csv');

maxIC50     = 1000;   % training peptides with a higher IC50 (uM) are left out of the fit, as in the original script
minLeafSize = [];     % empty = fitrtree default. The paper reports 570 (Bayesian optimisation); see README.

%% Train the regression tree
% ic50_train: 1466 x 1 IC50 values (uM); train_data_final: 1466 x 31 features
% (23 sequence features followed by 8 structure features).
load(trainFile, 'ic50_train', 'train_data_final');

[ic50_sorted, order] = sort(ic50_train);
Xtrain = train_data_final(order, :);
use    = ic50_sorted <= maxIC50;

rng default
if isempty(minLeafSize)
    model = fitrtree(Xtrain(use, :), ic50_sorted(use));
else
    model = fitrtree(Xtrain(use, :), ic50_sorted(use), 'MinLeafSize', minLeafSize);
end
fprintf('Trained on %d of %d peptides (IC50 <= %g uM).\n', nnz(use), numel(ic50_sorted), maxIC50);

%% Read the test peptides
% The two files list the same peptides in different orders, so they are matched by name.
[namesStruct, S] = read_table_csv(structFile,  true,  'first');   % name, then 8 features
[namesPseaac, P] = read_table_csv(pseaacFile, false, 'last');     % 23 features, then name

[found, loc] = ismember(namesPseaac, namesStruct);
if ~all(found)
    error('These peptides are in %s but not in %s: %s', pseaacFile, structFile, strjoin(namesPseaac(~found), ', '));
end

Xtest = [P, S(loc, :)];                      % same column order as the training matrix
if size(Xtest, 2) ~= size(Xtrain, 2)
    error('Test data has %d feature columns but the model was trained on %d.', size(Xtest, 2), size(Xtrain, 2));
end

[names, sortIdx] = sort(namesPseaac);
Xtest = Xtest(sortIdx, :);

%% Predict
ic50_predicted = abs(predict(model, Xtest));

% Classes as defined in the paper: positive <= 100 uM, negative > 1000 uM, otherwise average.
activityClass = repmat({'average'}, numel(ic50_predicted), 1);
activityClass(ic50_predicted <= 100)  = {'positive'};
activityClass(ic50_predicted > 1000)  = {'negative'};

result = table(names, ic50_predicted, activityClass, 'VariableNames', {'Peptide', 'IC50_predicted_uM', 'Class'});
disp(result)
writetable(result, outFile);
fprintf('Predictions written to %s\n', outFile);


function [names, X] = read_table_csv(file, hasHeader, nameColumn)
% Read a comma-separated file of peptide names and numeric features.
% nameColumn is 'first' or 'last'. Spaces around the fields are ignored.
txt   = strtrim(fileread(file));
lines = strsplit(txt, {'\r\n', '\n', '\r'});
if hasHeader
    lines(1) = [];
end
lines = lines(~cellfun(@(l) isempty(strtrim(l)), lines));

n     = numel(lines);
nCols = numel(strsplit(lines{1}, ','));
names = cell(n, 1);
X     = zeros(n, nCols - 1);
for i = 1:n
    p = strtrim(strsplit(lines{i}, ','));
    if numel(p) ~= nCols
        error('Line %d of %s has %d fields, expected %d.', i, file, numel(p), nCols);
    end
    if strcmp(nameColumn, 'first')
        names{i} = p{1};
        X(i, :)  = str2double(p(2:end));
    else
        names{i} = p{end};
        X(i, :)  = str2double(p(1:end-1));
    end
end
if any(isnan(X(:)))
    error('Non-numeric value found in %s.', file);
end
end
