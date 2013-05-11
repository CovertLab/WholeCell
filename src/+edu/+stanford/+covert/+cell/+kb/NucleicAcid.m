% Defines a nucleic acid. Base class for
% - dsDNA
%   - Genome
% - dsRNA
% - ssDNA
%   - Gene
%   - TranscriptionUnit
% - ssRNA
%   - mRNA
%   - rRNA
%   - tRNA
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/7/2009
classdef NucleicAcid < edu.stanford.covert.cell.kb.Polymer
    %computed properties
    properties %(SetAccess = protected)
        pI
        gcContent
    end

    methods
        function this = NucleicAcid(knowledgeBase, wid, wholeCellModelID, name,...
                sequence,...
                comments, crossReferences)

            if nargin == 0; return; end;

            this = edu.stanford.covert.cell.kb.NucleicAcid.empty(size(wid,1),0);
            this(size(wid,1),1) = edu.stanford.covert.cell.kb.NucleicAcid;
            for i=1:size(wid,1)
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
                        for j=1:size(fields,1)
                            values = crossReferences.(fields{j});
                            this(i,1).crossReferences.(fields{j}) = values(i);
                        end
                    end
                end

                this(i,1).sequence = sequence(i);
            end
        end

        %http://isoelectric.ovh.org/files/practise-isoelectric-point.html
        function value = get.pI(this)
            %retrieve
            if ~isempty(this.pI)
                value = this.pI;
                return;
            end
            
            %compute
            sequence = this.sequence;
            numA = sum('A' == sequence);
            numC = sum('C' == sequence);
            numG = sum('G' == sequence);
            numT = sum('T' == sequence);
            numU = sum('U' == sequence);

            pH = 6.5;             %starting point pI = 6.5 - theoretically it should be 7, but average protein pI is 6.5 so we increase the probability
            pHprev = 0.0;         %of finding the solution
            pHnext = 14.0;        %0-14 is possible pH range
            E = 0.01;             %epsilon means precision [pI = pH ± E]

            %the infinite loop
            while true

                % http://www.steve.gb.com/science/nucleic_acids.html
                NQ = 0;
                %NQ = NQ-1/(1+10^(3.65-pH));     %3' charge
                NQ = NQ-numA/(1+10^(3.5-pH));   %A charge
                NQ = NQ-numC/(1+10^(4.2-pH));   %C charge
                NQ = NQ+numG/(1+10^(pH-9.2));   %G charge
                NQ = NQ-numG/(1+10^(1.6-pH));   %G charge
                NQ = NQ+numT/(1+10^(pH-9.7));   %T charge
                NQ = NQ+numU/(1+10^(pH-9.2));   %U charge
                %NQ = NQ+1/(1+10^(pH-8.2));      %5' charge

                if pH >= 14.0
                    throw(MException('NucleicAcid:error','pH higher than 14'));
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%   BISECTION   %%%%%%%%%%%%%%%%%%%%%%%%

                %we are out of range, thus the new pH value must be smaller
                if NQ < 0
                    temp = pH;
                    pH = pH-((pH-pHprev)/2);
                    pHnext = temp;

                %we used to small pH value, so we have to increase it
                else
                    temp = pH;
                    pH = pH + ((pHnext-pH)/2);
                    pHprev = temp;
                end

                %terminal condition, finding isoelectric point with given precision
                if (pH-pHprev<E) && (pHnext-pH<E)
                    break;
                end
            end

            value = pH;
            
            %store
            this.pI = value;
        end

        function value = get.gcContent(this)
            %retrieve
            if ~isempty(this.gcContent)
                value = this.gcContent;
                return;
            end
            
            %compute
            sequence = this.sequence;
            value = (sum('G' == sequence) + sum('C' == sequence))/...
                this.sequenceLength;
            
            %store
            this.gcContent = value;
        end
    end
end