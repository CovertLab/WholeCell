% Run a simulation with modified enzyme kinetics
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/24/2013
classdef ModifiedMetabolicEnzymeKinetics < edu.stanford.covert.cell.sim.runners.SimulationRunner
    properties
        %Structure containing modified metabolic enzyme kinetics.
        %kCatFor, kCatRev should have units (rxn/enzyme/s).
        %Set kCatFor, kCatRev to NaN to use default value.
        modifiedKinetics = repmat(struct(...
            'reactionWholeCellModelID', '', ...
            'kCatFor', NaN, ...
            'kCatRev', NaN ...
            ), 0, 1); 
    end
    
    methods
        function this = ModifiedMetabolicEnzymeKinetics(varargin)
            this = this@edu.stanford.covert.cell.sim.runners.SimulationRunner(varargin{:});
        end
    end
    
    methods (Access = protected)
        function modifyNetworkParameters(this, sim)
            %get handles
            met = sim.process('Metabolism');
            
            %calculate indices of modified reactions
            ids = {this.modifiedKinetics.reactionWholeCellModelID}';
            kCatFor = [this.modifiedKinetics.kCatFor]';
            kCatRev = [this.modifiedKinetics.kCatRev]';
            
            if numel(unique(ids)) < numel(ids)
                throw(MException('ModifiedMetabolicEnzymeKinetics:invalidReactionList', ...
                    'Modified reaction list cannot contain multiple entries for a given reaction.'));
            end
            
            idxs = met.reactionIndexs(ids);
            
            %modify kcats
            met.enzymeBounds(idxs(~isnan(kCatRev)), 1) = kCatRev(~isnan(kCatRev));
            met.enzymeBounds(idxs(~isnan(kCatFor)), 2) = kCatFor(~isnan(kCatFor));
            met.fbaEnzymeBounds(met.fbaReactionIndexs_metabolicConversion, :) = met.enzymeBounds(met.reactionIndexs_fba, :);
        end
    end
end