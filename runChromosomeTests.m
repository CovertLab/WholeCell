% Run chromosome tests.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/13/2010

%% initialize
setWarnings();
setPath();
setPreferences();

%% run tests
runtests('-verbose', {
    'edu.stanford.covert.cell.sim.state.Chromosome_Test'
    'edu.stanford.covert.cell.sim.process.ChromosomeCondensation_Test'
    'edu.stanford.covert.cell.sim.process.ChromosomeSegregation_Test'
    'edu.stanford.covert.cell.sim.process.Cytokinesis_Test'
    'edu.stanford.covert.cell.sim.process.DNADamage_Test'
    'edu.stanford.covert.cell.sim.process.DNARepair_Test'
    'edu.stanford.covert.cell.sim.process.DNASupercoiling_Test'
    'edu.stanford.covert.cell.sim.process.Replication_Test'
    'edu.stanford.covert.cell.sim.process.ReplicationInitiation_Test'
    'edu.stanford.covert.cell.sim.process.Transcription_Test'
    'edu.stanford.covert.cell.sim.process.TranscriptionalRegulation_Test'
    'edu.stanford.covert.cell.sim.DNA_Test'
    'edu.stanford.covert.cell.sim.DNADamageRepair_Test'
    });