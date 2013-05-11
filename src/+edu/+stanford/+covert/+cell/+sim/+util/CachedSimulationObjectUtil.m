% Stores/loads simulatons
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 6/30/2011
classdef CachedSimulationObjectUtil
    methods (Static = true)
        function store(simulation, knowledgeBaseWID, fileName) %#ok<INUSL,INUSD>
            if nargin < 3
                fileName = 'data/Simulation_fitted%s.mat';
            end
            
            save(sprintf(fileName, ''), 'simulation', 'knowledgeBaseWID');
            
            revision = edu.stanford.covert.util.revision;
            if ~isnan(revision)
                save(sprintf(fileName, sprintf('-R%04d', revision + 1)), 'simulation', 'knowledgeBaseWID');
            end
        end
        
        function [sim, kbWID] = load(revision, fileName)
            if nargin < 2
                fileName = 'data/Simulation_fitted%s.mat';
            end
            
            if nargin == 0
                tmp = load(sprintf(fileName, ''));
            else
                files = dir(sprintf(fileName, '-R*'));
                fileRevisions = cellfun(@(name) str2double(name(20:23)), {files.name}');
                idx = find(fileRevisions <= revision, 1, 'last');
                if isempty(idx)
                    throw(MException('CachedSimulationObjectUtil:error', 'No cached simulation object matches revision %d', revision));
                end
                tmp = load(['data' filesep files(idx).name]);
            end
            sim = tmp.simulation;
            kbWID = tmp.knowledgeBaseWID;
            
            sim.constructRandStream();
        end
    end
end
