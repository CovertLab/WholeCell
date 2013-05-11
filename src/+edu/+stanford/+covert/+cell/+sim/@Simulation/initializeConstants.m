function this = initializeConstants(this, knowledgeBase)
% Cache various properties of knowledgeBase
% - mRNA, rRNA, tRNA gene indices within genes
% - ribosome 30S, 50S; RNA, DNA polymerase
%   within protein complexs
% - transcription unit properties
%   - binding probability
%   - sequence
%   - length
%   - gene composition
% - indices of transcription factors within protein monomers
%   - sigma factor
%   - elongation factor
%   - termination factor
% - indices of translation factors within protein monomers
%   - initiation factor
%   - elongation factors
%   - termination factors
% - gene properties
%   - RNA weight
%   - RNA base counts
%   - RNA half life
% - genetic code
%   - tRNA sequence
%   - tRNA aminoacylation
% - tRNA synthetases
%   - indices within protein monomers and protein complexs
%   - rates
% - composition of protein complexs
% - weights
%   - RNA
%   - protein monomers
%   - protein complexs
% - metabolism
%   - stoichiometry of reactions
%   - upper and lower transport bounds
%   - upper and lower enzyme bounds
%   - FBA right-hand side
%   - biomass composition
%   - FBA objective
%   - reaction catalysis matrix
%   - indices of dNTPs, NTPs, energy within metabolites
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/7/2011

%% constants
this.compartment.initializeConstants(knowledgeBase, this);
this.gene.initializeConstants(knowledgeBase, this);

%% states
this.state('Mass').initializeConstants(knowledgeBase, this);
this.state('Geometry').initializeConstants(knowledgeBase, this);
this.state('Stimulus').initializeConstants(knowledgeBase, this);
this.state('Metabolite').initializeConstants(knowledgeBase, this);
this.state('Rna').initializeConstants(knowledgeBase, this);
this.state('ProteinMonomer').initializeConstants(knowledgeBase, this);
this.state('ProteinComplex').initializeConstants(knowledgeBase, this);

states = this.states;
for i = 1:numel(states)
    states{i}.initializeConstants(knowledgeBase, this);
end

%% Processes
processes = this.processes;
for i = 1:length(processes)
    processes{i}.initializeConstants(knowledgeBase, this);
end
