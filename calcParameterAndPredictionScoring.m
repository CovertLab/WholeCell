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
% - submissionList
% - distances
% - Scores [double array]: score of each submission
% - Ranks [int array]: rank of each submission
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/12/2013
function [submissionList, distances, pValues, scores, ranks] = calcParameterAndPredictionScoring(varargin)
[submissionList, distances, pValues, scores, ranks]  = edu.stanford.covert.cell.sim.util.DreamScoring.calcScores(varargin{:});