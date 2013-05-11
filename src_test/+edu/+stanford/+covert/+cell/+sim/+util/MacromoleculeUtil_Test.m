% MacromoleculeUtility test
%
% Author: Derek Macklin, macklin@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated 6/29/2011
classdef MacromoleculeUtil_Test < TestCase
    properties
        simulation
    end
    
    methods
        function this = MacromoleculeUtil_Test(name)
            this = this@TestCase(name);
        end
        
        function setUp(this)
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load();
            this.simulation = sim;            
        end
        
        function testGetMacromoleculeIndexsIDsNames(this)
            import edu.stanford.covert.cell.sim.util.MacromoleculeUtil;
            
            s = this.simulation;
            
            % Gene that doesn't exist
            [idx, compIdx, id, name, type, ...
                geneIdxs, geneCompIdxs, geneIds, geneNames, ...
                nascentRnaIdxs, nascentRnaCompIdxs, nascentRnaIds, nascentRnaNames, ...
                matureRnaIdxs, matureRnaCompIdxs, matureRnaIds, matureRnaNames, ...
                monomerIdxs, monomerCompIdxs, monomerIds, monomerNames, ...
                complexIdxs, complexCompIdxs, complexIds, complexNames] = ...
                MacromoleculeUtil.getMacromoleculeIndexsIDsNames('DOES_NOT_EXIST', s);
            
            assert(isequal(idx,[]));
            assert(isequal(compIdx,[]));
            assert(isequal(id,[]));
            assert(isequal(name,[]));
            assert(isequal(type,[]));
            assert(isequal(geneIdxs,[]));
            assert(isequal(geneCompIdxs,[]));
            assert(isequal(geneIds,[]));
            assert(isequal(geneNames,[]));
            assert(isequal(nascentRnaIdxs,[]));
            assert(isequal(nascentRnaCompIdxs,[]));
            assert(isequal(nascentRnaIds,[]));
            assert(isequal(nascentRnaNames,[]));
            assert(isequal(matureRnaIdxs,[]));
            assert(isequal(matureRnaCompIdxs,[]));
            assert(isequal(matureRnaIds,[]));
            assert(isequal(matureRnaNames,[]));
            assert(isequal(monomerIdxs,[]));
            assert(isequal(monomerCompIdxs,[]));
            assert(isequal(monomerIds,[]));
            assert(isequal(monomerNames,[]));
            assert(isequal(complexIdxs,[]));
            assert(isequal(complexCompIdxs,[]));
            assert(isequal(complexIds,[]));
            assert(isequal(complexNames,[]));
            
            % Two species present in a protein complex
            [idx, compIdx, id, name, type, ...
                geneIdxs, geneCompIdxs, geneIds, geneNames, ...
                nascentRnaIdxs, nascentRnaCompIdxs, nascentRnaIds, nascentRnaNames, ...
                matureRnaIdxs, matureRnaCompIdxs, matureRnaIds, matureRnaNames, ...
                monomerIdxs, monomerCompIdxs, monomerIds, monomerNames, ...
                complexIdxs, complexCompIdxs, complexIds, complexNames] = ...
                MacromoleculeUtil.getMacromoleculeIndexsIDsNames('MG_008_379_TETRAMER',s);
            
            assert(isequal(id,'MG_008_379_TETRAMER'));
            
            assert(isequal(name,'tRNA uridine 5-carboxymethylaminomethyl modification enzyme'));
            
            assert(isequal(type,'complex'));
            
            assert(isequal(size(geneIdxs),[2 1]));
            
            assert(isequal(size(geneCompIdxs),[2 1]));
            
            assert(isequal(size(geneIds),[2 1]));
            assert(isequal(geneIds,{'MG_008';'MG_379'}));
            
            assert(isequal(size(geneNames),[2 1]));
            assert(isequal(geneNames,{'tRNA modification GTPase TrmE';'tRNA uridine 5-carboxymethylaminomethyl modification enzyme GidA'}));
            
            
            assert(isequal(size(nascentRnaIdxs),[1 2]));
            
            assert(isequal(size(nascentRnaCompIdxs),[1 2]));
            
            assert(isequal(size(nascentRnaIds),[2 1]));
            assert(isequal(nascentRnaIds,{'TU_003';'TU_275'}));
            
            
            assert(isequal(size(nascentRnaNames),[2 1]));
            assert(isequal(nascentRnaNames,{'MG_005 - MG_009';'MG_379 - MG_381'}));
            
            
            assert(isequal(size(matureRnaIdxs),[1 2]));
            
            assert(isequal(size(matureRnaCompIdxs),[1 2]));
            
            assert(isequal(size(matureRnaIds),[2 1]));
            assert(isequal(matureRnaIds,{'TU_003';'TU_275'}));
            
            assert(isequal(size(matureRnaNames),[2 1]));
            assert(isequal(matureRnaNames,{'MG_005 - MG_009';'MG_379 - MG_381'}));
            
            assert(isequal(size(monomerIdxs),[2 1]));
            assert(isequal(size(monomerCompIdxs),[2 1]));
            
            assert(isequal(size(monomerIds),[2 1]));
            assert(isequal(monomerIds,{'MG_008_MONOMER';'MG_379_MONOMER'}));
            
            assert(isequal(size(monomerNames),[2 1]));
            assert(isequal(monomerNames,{'tRNA modification GTPase TrmE';'tRNA uridine 5-carboxymethylaminomethyl modification enzyme GidA'}));
            
            assert(isequal(complexIds,{'MG_008_379_TETRAMER'}));
            assert(isequal(complexNames,{'tRNA uridine 5-carboxymethylaminomethyl modification enzyme'}));
        end
        
        function testGetMacromoleculeCounts(~)            
            import edu.stanford.covert.cell.sim.util.MacromoleculeUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            %select a simulation to analyze
            [simDir, ~, sim] = SimulationDiskUtil.getLatestSimulation();
            
            %load data
            [cnts, ...
                geneCnts, ...
                nascentRnaCnts, processedRnaCnts, matureRnaCnts, boundRnaCnts, misfoldedRnaCnts, damagedRnaCnts, aminoacylatedRnaCnts, ...
                nascentMonCnts, processedIMonCnts, translocatedMonCnts, processedIIMonCnts, foldedMonCnts, matureMonCnts, inactivatedMonCnts, boundMonCnts, misfoldedMonCnts, damagedMonCnts, ...
                nascentCpxCnts, matureCpxCnts, inactivatedCpxCnts, boundCpxCnts, misfoldedCpxCnts, damagedCpxCnts, ...
                rnaPolTUs, rnaPolPosStrnds, ribMRNAs, ribPos] = ...
                MacromoleculeUtil.getMacromoleculeCounts('MG_357', sim, simDir);
        end
    end
end

