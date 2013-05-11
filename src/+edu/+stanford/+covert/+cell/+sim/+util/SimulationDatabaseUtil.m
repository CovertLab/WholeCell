%SimulationDatabaseUtil
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 6/30/2011
classdef SimulationDatabaseUtil
    methods (Static = true)
        function WID = selectSimulation(knowledgeBaseWID, database)
            database.setNullValue(0);
            database.prepareStatement('CALL get_simulations("{Si}")', knowledgeBaseWID);
            simulations = database.query();            
            
            %display available simulations
            fprintf('%s\t%s\t\t%s\t%s\t\t%s\n','#','Description             ','Revision','Date               ','Investigator');
            fprintf('%s\t%s\t\t%s\t%s\t\t%s\n','=','========================','========','===================','============');
            for i=1:length(simulations.WID)
                shortDescription=[simulations.ShortDescription{i} '             '];
                shortDescription=shortDescription(1:min(31,length(shortDescription)));
                fprintf(['%d.\t%s' repmat('\t',1,floor((37-length(shortDescription))/4)) '%d\t\t\t%s\t\t%s\n'],...
                    i,shortDescription,simulations.Revision(i),simulations.StartDate{i},simulations.Investigator{i});
            end
            
            %select knowledge base
            selectedSimulation = 0;
            while ~isa(selectedSimulation,'double') || mod(selectedSimulation,1) || ...
                    selectedSimulation < 1 || selectedSimulation > length(simulations.WID)
                selectedSimulation = input('Please select a simulation (by number): ');
            end
            
            %return WID
            WID = simulations.WID(selectedSimulation);
        end
    end
end