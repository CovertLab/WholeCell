% Defines a dsDNA polymer. Base class for
% - Genome
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/7/2009
classdef dsDNA < edu.stanford.covert.cell.kb.NucleicAcid
    properties %(SetAccess = protected)
        empiricalFormula
        smiles
        charge
        pKa
    end
    
    %computed properties
    properties %(SetAccess = protected)
        baseCount
        cumulativeBaseCount
        decayReaction
        molecularWeight
        density
        volume
        extinctionCoefficient
        absorbanceFactor
    end
    
    properties (Constant = true)
        singleExtinction = [15400 7400 11500 8700];
        pairwiseExtinction = [
            27400 21200 25000 22800;
            21200 14600 18000 15200;
            25200 17600 21600 20000;
            23400 16200 19000 16800];
    end
    
    methods
        function this = dsDNA(knowledgeBase, wid, wholeCellModelID, name, ...
                sequence, ...
                comments, crossReferences)
            
            if nargin == 0; return; end;
            
            this = edu.stanford.covert.cell.kb.dsDNA.empty(size(wid,1),0);
            this(size(wid,1),1) = edu.stanford.covert.cell.kb.dsDNA;
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
                
                this(i,1).sequence = sequence(i);
            end
        end
        
        function value = get.empiricalFormula(~)
            throw(MException('dsDNA:error', 'property is not defined'));
        end
        
        function value = get.smiles(~)
            throw(MException('dsDNA:error', 'property is not defined'));
        end
        
        function value = get.charge(~)
            throw(MException('dsDNA:error', 'property is not defined'));
        end
        
        function value = get.pKa(~)
            throw(MException('dsDNA:error', 'property is not defined'));
        end
        
        function value = get.baseCount(this)
            %retrieve
            if ~isempty(this.baseCount)
                value = this.baseCount;
                return;
            end
            
            %calculate
            value = this.computeBaseCount(this.sequence, this.knowledgeBase.numMetabolites, this.knowledgeBase.dnmpIndexs);
            
            %store
            this.baseCount = value;
        end
        
        function value = get.cumulativeBaseCount(this)
            %retrieve
            if ~isempty(this.cumulativeBaseCount)
                value = this.cumulativeBaseCount;
                return;
            end
            
            %calculate
            sequence = this.sequence';
            value = zeros(length(sequence),this.knowledgeBase.numMetabolites);
            value(:,this.knowledgeBase.dnmpIndexs) = 2*[...
                sequence == 'A' ...
                sequence == 'C' ...
                sequence == 'G' ...
                sequence == 'T'];
            
            value = cumsum(value,2);
            
            %store
            this.cumulativeBaseCount = value;
        end
        
        function value = get.decayReaction(this)
            %retrieve
            if ~isempty(this.decayReaction)
                value = this.decayReaction;
                return;
            end
            
            %calculate
            value = this.computeDecayReaction(this.baseCount, this.sequenceLength, this.sequenceTopology, this.knowledgeBase.waterIndexs, this.knowledgeBase.hydrogenIndexs);
            
            %store
            this.decayReaction = value;
        end
        
        function value = get.molecularWeight(this)
            %retrieve
            if ~isempty(this.molecularWeight)
                value = this.molecularWeight;
                return;
            end
            
            %calculate
            value = this.calculateMolecularWeight(...
                this.sequence, this.sequenceLength, this.sequenceTopology, ...
                this.knowledgeBase.metaboliteMolecularWeights(this.knowledgeBase.dnmpIndexs));
            
            %store
            this.molecularWeight = value;
        end
        
        function value = get.density(this)
            %retrieve
            if ~isempty(this.density)
                value = this.density;
                return;
            end
            
            %calculate
            value = this.molecularWeight/this.volume;
            
            %store
            this.density = value;
        end
        
        %http://www.ncbi.nlm.nih.gov/pubmed/6639645
        function value = get.volume(this)
            %retrieve
            if ~isempty(this.volume)
                value = this.volume;
                return;
            end
            
            %calculate
            value = pi * (12e-7)^2 * 0.34e-7 * this.sequenceLength*edu.stanford.covert.util.ConstantUtil.nAvogadro;
            
            %store
            this.volume = value;
        end
        
        %Extinction (absorption) coefficients at 260nm
        %http://www.owczarzy.net/extinct.htm
        function value = get.extinctionCoefficient(this)
            %retrieve
            if ~isempty(this.extinctionCoefficient)
                value = this.extinctionCoefficient;
                return;
            end
            
            %calculate
            sequence = this.sequence;
            idx2 = nt2int(sequence(1), 'Alphabet', 'DNA');
            value = 0;
            for i = 1:length(sequence)-1
                idx1 = idx2;
                idx2 = n2int(sequence(2), 'Alphabet', 'DNA');
                if idx1 > 4 || idx2 > 4; continue; end;
                value = value + this.pairwiseExtinction(idx1, idx2) - this.singleExtinction(idx2);
            end
            value = value + this.singleExtinction(idx2);
            
            %store
            this.extinctionCoefficient = value;
        end
        
        %absorbance factor (mmol^-1)
        function value = get.absorbanceFactor(this)
            %retrieve
            if ~isempty(this.absorbanceFactor)
                value = this.absorbanceFactor;
                return;
            end
            
            %calculate
            value = 1 / this.extinctionCoefficient;
            
            %store
            this.absorbanceFactor = value;
        end
    end
    
    methods (Static = true)
        function value = computeBaseCount(sequence, numMetabolites, dnmpIndexs)
            value = zeros(1, numMetabolites);
            value(dnmpIndexs) = 2 * [...
                sum('A' == sequence) ...
                sum('C' == sequence) ...
                sum('G' == sequence) ...
                sum('T' == sequence)];
        end
        
        function value = computeDecayReaction(baseCount, sequenceLength, sequenceTopology, waterIndexs, hydrogenIndexs)
            value = baseCount;
            value(waterIndexs)    = value(waterIndexs)    - 2 * max(0, sequenceLength - strcmp(sequenceTopology, 'linear'));
            value(hydrogenIndexs) = value(hydrogenIndexs) + 2 * max(0, sequenceLength - strcmp(sequenceTopology, 'linear'));
        end
        
        function value = calculateMolecularWeight(sequence, sequenceLength, sequenceTopology, dnmpMolecularWeights)
            % import classes
            import edu.stanford.covert.util.ConstantUtil;
                      
            molecularWeightHO = ConstantUtil.elements.H + ConstantUtil.elements.O;
            
            dnmpCount = [
                nnz(sequence == 'A');
                nnz(sequence == 'C');
                nnz(sequence == 'G');
                nnz(sequence == 'T')];
            dnmpCount = dnmpCount + dnmpCount([4 3 2 1]);
            
            value = dnmpCount' * dnmpMolecularWeights ...
                - 2 * molecularWeightHO * (sequenceLength - strcmp(sequenceTopology, 'linear'));
        end
    end
end