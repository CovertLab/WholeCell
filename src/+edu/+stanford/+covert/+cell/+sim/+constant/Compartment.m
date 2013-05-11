classdef Compartment < handle
    properties
        count             %number of compartments
        wholeCellModelIDs %ids of compartments
        names             %names of compartments
    end
    
    %indices
    properties
        cytosolIndexs                   %index within compartments
        chromosomeIndexs                %index within compartments
        membraneIndexs                  %index within compartments
        terminalOrganelleCytosolIndexs  %index within compartments
        terminalOrganelleMembraneIndexs %index within compartments
        extracellularIndexs             %index within compartments
        cellularIndexs                  %index within compartments
    end
    
    methods
        function this = Compartment()
        end
        
        function initializeConstants(this, knowledgeBase, ~)
            this.count             = knowledgeBase.numCompartments;
            this.wholeCellModelIDs = {knowledgeBase.compartments.wholeCellModelID}';
            this.names             = {knowledgeBase.compartments.name}';
            
            this.cytosolIndexs                   = knowledgeBase.cytosolCompartmentIndexs;
            this.chromosomeIndexs                = knowledgeBase.chromosomeCompartmentIndexs;
            this.membraneIndexs                  = knowledgeBase.membraneCompartmentIndexs;
            this.terminalOrganelleCytosolIndexs  = knowledgeBase.terminalOrganelleCytosolCompartmentIndexs;
            this.terminalOrganelleMembraneIndexs = knowledgeBase.terminalOrganelleMembraneCompartmentIndexs;
            this.extracellularIndexs             = knowledgeBase.extracellularCompartmentIndexs;
            this.cellularIndexs                  = knowledgeBase.cellularCompartmentIndexs;            
        end
        
        function index = getIndexs(this, wholeCellModelIDs)
            [~, index] = ismember(wholeCellModelIDs, this.wholeCellModelIDs);
        end
    end
end