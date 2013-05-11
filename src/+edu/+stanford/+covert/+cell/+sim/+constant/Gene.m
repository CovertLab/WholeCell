classdef Gene < handle
    %indices of cellular constituents
    properties (SetAccess = private)
        wholeCellModelIDs     %Whole Cell model IDs
        names                 %names
        startCoordinates      %start coordinates
        lengths               %lengths
        strands               %strands
        essential             %is gene essential (Y=>yes, N=>no, M=>maybe)
        
        mRNAIndexs            %mRNA genes within geneWholeCellModelIDs
        rRNAIndexs            %rRNA genes within geneWholeCellModelIDs
        sRNAIndexs            %sRNA genes within geneWholeCellModelIDs
        tRNAIndexs            %tRNA genes within geneWholeCellModelIDs
        ribosomalRRNAIndexs   %ribosomal RNA genes within geneWholeCellModelIDs
    end
    
    methods
        function this = Gene()
        end
        
        function initializeConstants(this, knowledgeBase, ~)
            this.wholeCellModelIDs   = {knowledgeBase.genes.wholeCellModelID}';
            this.names               = {knowledgeBase.genes.name}';
            this.startCoordinates    = [knowledgeBase.genes.startCoordinate]';
            this.lengths             = [knowledgeBase.genes.sequenceLength]';
            this.strands             = 2-[knowledgeBase.genes.direction]';
            this.essential           = {knowledgeBase.genes.essential}';
            
            this.mRNAIndexs          = double([knowledgeBase.mRNAGenes.idx])';
            this.rRNAIndexs          = double([knowledgeBase.rRNAGenes.idx])';
            this.sRNAIndexs          = double([knowledgeBase.sRNAGenes.idx])';
            this.tRNAIndexs          = double([knowledgeBase.tRNAGenes.idx])';
            this.ribosomalRRNAIndexs = double(knowledgeBase.ribosomalRRNAIndexs);
        end
    end
    
    %helper functions
    methods
        function value = getIndexs(this, wholeCellModelIDs)
            [~, value] = ismember(wholeCellModelIDs, this.wholeCellModelIDs);
        end
        
        function value = getIndexsByPosition(this, pos)
            value = find(...
                pos >= this.startCoordinates & ...
                pos <= this.startCoordinates + this.lengths - 1 ...
                , 1);
        end
    end
end