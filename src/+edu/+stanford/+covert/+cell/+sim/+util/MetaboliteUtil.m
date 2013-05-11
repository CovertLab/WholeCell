classdef MetaboliteUtil
    methods (Static)
        function localizedMetabolites = localizeMetabolites(sim, metabolites)
            mets = sim.state('Metabolite');
            nMets = numel(mets.wholeCellModelIDs);
            nComp = sim.compartment.count;
            cIdx = sim.compartment.cytosolIndexs;
            mIdx = sim.compartment.membraneIndexs;
            hMetIdxs = mets.hydrophobicIndexs;
            
            localizedMetabolites = zeros(nMets, nComp);
            localizedMetabolites(:, cIdx) = metabolites;
            localizedMetabolites(hMetIdxs, mIdx) = localizedMetabolites(hMetIdxs, cIdx);
            localizedMetabolites(hMetIdxs, cIdx) = 0;
        end
    end
end