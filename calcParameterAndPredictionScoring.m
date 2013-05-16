%calcParameterAndPredictionScoring
% Calculates score and rank of each submission by the procedure established
% in the DREAM 7 parameter estimation challenge (see
% http://www.the-dream-project.org/sites/the-dream-project.org/files/documents/DREAM6/DREAM7/Descriptions/dream7c1_scoring.pdf)
%
% Input:
% - dists [struct array]: struct array of distances calculated by
%   calcParametersAndPredictionError
%
% Output
% - Scores [double array]: score of each submission
% - Ranks [int array]: rank of each submission
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/12/2013

% TODO: get dists from Synapse
function [scores, ranks] = calcParameterAndPredictionScoring(dists)

p_param = sampleAndCalcPVal([dists.parameter]', 1000);
p_pred = sampleAndCalcPVal([dists.prediction]', 1000);

scores = -log(p_param .* p_pred);

[~, idxs] = sort(scores);
ranks = zeros(size(scores));
ranks(idxs) = numel(scores):-1:1;

function p = sampleAndCalcPVal(dists, nTrial)
nDists = numel(dists);
p = zeros(nDists, nTrial);
for iTrial = 1:nTrial
    distrib = dists(ceil(rand(nDists, 1) * nDists));
    p(:, iTrial) = calcPVal(dists, distrib);    
end

%average over trials
p = mean(p, 2);


%calculate p-value
function p = calcPVal(tests, distrib)
[sortedTestAndDistrib, idx] = sortrows([
    tests zeros(size(tests))
    distrib ones(size(distrib))
    ], [-1 -2]);
percentiles = zeros(size(sortedTestAndDistrib, 1), 1);
percentiles(idx) = cumsum(sortedTestAndDistrib(:, 2)) / numel(distrib);

p = percentiles(1:numel(tests));
p(p == 0) = 1 / numel(distrib);