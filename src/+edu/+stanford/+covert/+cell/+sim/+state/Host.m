%Host
%
% @wholeCellModelID State_Host
% @name             Host
% @description
%   This class reports qualitative metrics regarding the ability of M.
%   genitalium to interact with host human urogenital tract epithelial
%   cells:
%   - ability of M. genitalium to adhere to the urogenital epithelium,
%     which requires adhesins and a functional terminal organelle
%   - ability of M. genitalium lipoproteins to interact with host TLR
%     receptors 1, 2, and 6
%
%   References
%   ==========
%   See HostInteraction process.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/25/2011
classdef Host < edu.stanford.covert.cell.sim.CellState
    %property annotations
    properties (Constant)
        optionNames             = {   %names of properties that are options
            'verbosity';
            'seed';
            };
        fixedConstantNames      = {}; %names of process properties that are considered fixed constants
        fittedConstantNames     = {}; %names of process properties that are considered fitted constants, and should be stored with the simulation as such
        stateNames              = {   %names of properties which are part of the simulation's state
            'isBacteriumAdherent'
            'isTLRActivated'
            'isNFkBActivated'
            'isInflammatoryResponseActivated'
            };
        dependentStateNames     = {}; %names of properties which can be calculated from the simulation's state
    end
    
    properties (Constant = true)
        nTLRs = 3;
        tlrIDs = {
            'TLR1'
            'TLR2'
            'TLR6'
            };
        tlrIndexs_1 = 1;
        tlrIndexs_2 = 2;
        tlrIndexs_6 = 3;
    end
    
    %state
    properties
        isBacteriumAdherent             %true/false indicating whether or not cell can adhere
        isTLRActivated                  %true/false indicating the state of TLRs 1, 2, and 6 activated
        isNFkBActivated                 %true/false indicating the state of NF-kB
        isInflammatoryResponseActivated %true/false indicating the state of the inflammatory response
    end
    
    properties (Constant)
        dryWeight = 0; %dry weight of this class' state properties
    end
    
    %constructor
    methods
        function this = Host(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.CellState(wholeCellModelID, name);
        end
    end
    
    %allocate memory for state
    methods
        function allocateMemory(this, numTimePoints)
            this.isBacteriumAdherent = false(1, 1, numTimePoints);
            this.isTLRActivated = false(this.nTLRs, 1, numTimePoints);
        end
    end
    
    %initialization
    methods
        %state initialized by host interaction process
        function initialize(this)
        end
    end
end
