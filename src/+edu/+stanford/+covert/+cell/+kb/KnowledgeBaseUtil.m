% Utility functions for working with the knowledge base database.
% - get WID of latest knowledge base, simulation
% - select knowledge base and simulation from lists of available knowledge
%   bases and simulations, and return WID
% - get archived code of a simulation
% - create, retrieve contacts
% - create connection to knowledge base
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/7/2009
classdef KnowledgeBaseUtil
    methods (Static = true)
        %Get WID of the latest knowledge base
        function WID = selectLatestKnowledgeBase(database)
            database.setNullValue(0);
            database.prepareStatement('CALL get_latest_knowledgebase(''WholeCell'')');
            WID = database.query().WID(1);
        end

        %Prompt user to select a knowledge base from the list of avaiable
        %knowledge bases
        function WID = selectKnowledgeBase(database)
            %query available datasets
            database.setNullValue(0);
            database.prepareStatement('CALL get_knowledgebases("{S}","{S}")', 'WholeCell', 'Y');
            knowledgebases = database.query();

            %display available datasets
            fprintf('%s\t%s\t\t\t%s\t\t%s\t\t%s\n', '#', 'Name', 'Version', '    Insert Date    ', 'Investigator');
            fprintf('%s\t%s\t\t\t%s\t\t%s\t\t%s\n', '=', '====', '=======', '===================', '============');
            for i=1:length(knowledgebases.WID)
                fprintf(['%d.\t%s' repmat('\t',1,floor((19-length(knowledgebases.Name{i}))/4)) '%s\t\t%s\t\t%s\n'], ...
                    i, knowledgebases.Name{i}, knowledgebases.Version{i}, knowledgebases.InsertDate{i}, knowledgebases.Investigator{i});
            end

            %select knowledge base
            selectedKnowledgeBase = 0;
            while ~isa(selectedKnowledgeBase, 'double') || mod(selectedKnowledgeBase, 1) ~= 0 || ...
                    selectedKnowledgeBase < 1 || selectedKnowledgeBase > length(knowledgebases.WID)
                selectedKnowledgeBase = input('Please select a dataset (by number): ');
            end

            %return WID
            WID = knowledgebases.WID(selectedKnowledgeBase);
        end

        %Get WID of the latest simulation for knowledge base
        function simulationWID = selectLatestSimulation(knowledgeBaseWID, database)
            database.setNullValue(0);
            database.prepareStatement('CALL get_latest_simulation("{Si}")', knowledgeBaseWID);
            simulationWID = database.query().WID;
        end

        %Prompt user to select a simulation from the list of avaiable
        %simulations for the selected knowledge base
        function simulationWID = selectSimulation(knowledgeBaseWID, database)
            %retrieve list of simulations from database
            database.setNullValue(0);
            database.prepareStatement('CALL get_simulations("{Si}")', knowledgeBaseWID);
            simulations = database.query();

            %display list of available simulations
            fprintf('%s\t%s\t\t\t\t\t\t\t\t\t\t\t%s\t%s\t\t\t\t\t%s\n', '#', 'Label', 'Investigator', 'Date', 'Description');
            fprintf('%s\t%s\t\t\t\t\t\t\t\t\t\t\t%s\t%s\t\t\t\t\t%s\n', '=', '=====', '============', '====', '===========');
            for i = 1:length(simulations.WID)
                fprintf(['%d.\t%s' repmat('\t',1,ceil((48-length(simulations.Label{i}))/4)) '%s\t%s\t\t%s\n'], ...
                    i, simulations.Label{i}, simulations.Investigator{i}, simulations.Date{i}, char(simulations.Description{i}));
            end

            %prompt user to select simulation
            selectedSimulation=0;
            while ~isa(selectedSimulation, 'double') || mod(selectedSimulation, 1) ~= 0 || ...
                    selectedSimulation < 1 || selectedSimulation > length(simulations.WID)
                selectedSimulation = input('Please select a simulation: ');
            end

            %return WID
            simulationWID = simulations.WID(selectedSimulation);
        end

        %Retrieves zip archive of simulation stored in BioWhareouse, stores
        %zip archive to fileName
        function getSimulationCodeArchive(simulationWID, fileName, database)
            %retrieve archive from database
            database.setNullValue(0);
            database.prepareStatement('CALL get_simulation_codearchive("{Si}")', simulationWID);
            archive = database.query();

            %store archive to fileName
            fid = fopen(fileName, 'w+');
            fwrite(fid, cell2mat(archive.DifferencesFromRevision));
            fclose(fid);
        end

        function wid = getContact(investigator, knowledgeBaseWID, database)
            %get WID of contact
            database.setNullValue(0);
            database.prepareStatement('CALL get_contact("{S}","{Si}")', ...
                investigator.email, knowledgeBaseWID);
            wid = database.query().WID;

            %is no contact exists, create contact
            if isempty(wid)
                wid = edu.stanford.covert.cell.kb.KnowledgeBaseUtil.newContact(investigator, knowledgeBaseWID, database);
            end
        end

        %create contact
        function wid = newContact(investigator, knowledgeBaseWID, database)
            database.setNullValue(0);
            database.prepareStatement('CALL set_contact("{S}","{S}","{S}","{S}","{Si}")', ...
                investigator.firstName, investigator.lastName, investigator.email, investigator.affiliation, ...
                knowledgeBaseWID);
            wid = database.query().WID;
        end
    end
end