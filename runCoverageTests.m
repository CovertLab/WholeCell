% Gathers code coverage for whole cell simulation tests.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 12/15/2010
function runCoverageTests()
%% profile on
profile('on', '-nohistory');

%% run tests
setWarnings;
setPath();
setPreferences();
runSmallTests();
runMediumTests();
runMediumProteinTests();
runMediumDNATests();
runAnalysisTests();

%% profile off, save data
profile('off');
profData = profile('info'); %#ok<NASGU>
save('output/runCoverageTests/profData.mat', 'profData');

%% generate coverage report
report = edu.stanford.covert.test.Coverage('src', '..');
report.exportXML('output/runCoverageTests/coverage.xml');

%% generate documentation
edu.stanford.covert.cell.sim.util.Documentation.generate();
