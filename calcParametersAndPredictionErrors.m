%calcParametersAndPredictionErrors
% 
% Input
% - parameterVals [struct]: struct containing desired values of simulation
%   parameters. Initialize using the getAllParameters method of the
%   simulation class to get default parameter values. Overwrite struct
%   elements to set parameter values.
% - parameterValsPath [.mat/.xml file path]: .mat file path for stored
%   struct of parameter values or .xml file describing parameter values.
%   Use http://wholecell.stanford.edu/simulation/runSimulations.php to
%   generate XML file.
% - avgVals [struct]: Struct containing average value of in silico
%   experimental data
% - avgValsPath [.mat file path]: .mat file path to struct containing
%   average value of in silico experimental data
% - refParameterVals: same as parameterVals, but for reference parameter
%   values
% - refParameterValsPath: same as parameterValsPath, but for reference parameter
%   values
% - refAvgVals: same as avgVals, but for reference parameter
%   values
% - refAvgValsPath: same as avgValsPath:, but for reference parameter
%   values
%
% Output
% - dists [struct]: struct with two fields (parameter, prediction)
%   containing the cacluated parameter and prediction error from comparison
%   to the reference 
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/12/2013
function dists = calcParametersAndPredictionErrors(varargin)
dists = edu.stanford.covert.cell.sim.util.DreamScoring.calcErrors(varargin{:});