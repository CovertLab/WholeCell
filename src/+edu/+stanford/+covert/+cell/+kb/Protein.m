% Defines a protein. Base class for
% - ProteinMonomer
% - ProteinComplex
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/7/2009
classdef Protein < edu.stanford.covert.cell.kb.Enzyme
    %computed properties
    properties (Abstract = true) %SetAccess = protected
        negAA
        posAA

        instabilityIndex
        stable
        aliphaticIndex
        gravy               %GRAVY (Grand Average of Hydropathy)

        baseCount
        cumulativeBaseCount
        decayReaction
        tRNACount
    end
    
    properties (SetAccess = protected)
        dnaFootprint
        dnaFootprintBindingStrandedness
        dnaFootprintRegionStrandedness
    end

    methods
        function this = Protein(knowledgeBase, wid, wholeCellModelID, name,...
                molecularInteraction, chemicalRegulation, subsystem, ...
                generalClassification, proteaseClassification, ...
                transporterClassification, ...
                comments, crossReferences)

            if nargin == 0; return; end;

            this = edu.stanford.covert.cell.kb.Protein.empty(size(wid,1),0);
            this(size(wid,1),1) = edu.stanford.covert.cell.kb.Protein;
            for i = 1:size(wid,1)
                this(i,1).idx = i;
                this(i,1).knowledgeBase = knowledgeBase;
                this(i,1).wid = wid(i);
                this(i,1).wholeCellModelID = wholeCellModelID{i};
                this(i,1).name = name{i};
                if exist('comments','var') && ~isempty(comments); this(i,1).comments = comments{i}; end;
                if exist('crossReferences','var')
                    if size(crossReferences,1)>1
                        this(i,1).crossReferences = crossReferences(i);
                    else
                        this(i,1).crossReferences = struct;
                        fields = fieldnames(crossReferences);
                        for j = 1:size(fields,1)
                            values = crossReferences.(fields{j});
                            this(i,1).crossReferences.(fields{j}) = values(i);
                        end
                    end
                end

                this(i,1).molecularInteraction = molecularInteraction{i};
                this(i,1).chemicalRegulation = chemicalRegulation{i};
                this(i,1).subsystem = subsystem{i};
                this(i,1).generalClassification = generalClassification{i};
                this(i,1).proteaseClassification = proteaseClassification;
                this(i,1).transporterClassification = transporterClassification{i};
            end
        end
        
        %If DNA footprint hasn't been measured experimentally, estimates DNA
        %footprint from sequence:
        %- calculates molecular weight
        %- assumes density of 1.35 g/mol to calculate volume
        %- assumes protein is is sphere to calculate diameter
        %- assumes average length of base of DNA is 3.4 Angstroms to calculate
        %  the number of bases of DNA a protein spans
        %
        %units
        %- volume: ml/mol = 10^24 A^3/mol = 10^24/nAvogadro A^3 / molecule
        %- d: Angstroms
        function value = get.dnaFootprint(this)
            import edu.stanford.covert.util.ConstantUtil;
            
            if isempty(this.dnaFootprint) || isnan(this.dnaFootprint) || this.dnaFootprint < 1
                d = 2 * (3 / 4 / pi * this.volume / ConstantUtil.nAvogadro) ^ (1/3) * 1e8; 
                this.dnaFootprint = max(1, round(d / 3.4));
            end
            
            value = this.dnaFootprint;
        end
    end
end