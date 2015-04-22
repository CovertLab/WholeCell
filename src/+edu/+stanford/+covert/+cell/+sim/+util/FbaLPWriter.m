%FbaLPWriter
% Formats metabolic model in stm format.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updaetd 2/3/2011
classdef FbaLPWriter
    methods (Static)
        function write(sim, dirName)
            import edu.stanford.covert.cell.sim.util.FbaLPWriter;
            
            p = sim.process('Metabolism');
            
            %SBML-COBRA
            p.formulateFBA([], [], false);
            FbaLPWriter.write_cobra_matlab(sim, [dirName filesep 'metabolism.full.mmol-gDCW-h.sbml'], 'mmol/gDCW/h');
            FbaLPWriter.write_sbml(sim, [dirName filesep 'metabolism.full.mmol-gDCW-h.sbml'], 'mmol/gDCW/h');
            
            %LP
            p.formulateFBA([], [], false);
            FbaLPWriter.write_lpsolve_lp(sim, [dirName filesep 'metabolism.full.mmol-gDCW-h.lp'], 'mmol/gDCW/h');
            FbaLPWriter.write_cplex_lp(sim, [dirName filesep 'metabolism.full.mmol-gDCW-h.lpt'], 'mmol/gDCW/h');
            FbaLPWriter.write_gurobi_lp(sim, [dirName filesep 'metabolism.full.mmol-gDCW-h.gurobi.lp'], 'mmol/gDCW/h');
            
            %lindo
            p.formulateFBA([], [], false);
            FbaLPWriter.write_lindo(sim, [dirName filesep 'metabolism.full.molecules-cell-s.ltx'], 'molecules/cell/s');
            FbaLPWriter.write_lindo(sim, [dirName filesep 'metabolism.full.mmol-gDCW-h.ltx'], 'mmol/gDCW/h');
            
            p.formulateFBA([], [], true);
            FbaLPWriter.write_lindo(sim, [dirName filesep 'metabolism.reduced.molecules-cell-s.ltx'], 'molecules/cell/s');
            FbaLPWriter.write_lindo(sim, [dirName filesep 'metabolism.reduced.mmol-gDCW-h.ltx'], 'mmol/gDCW/h');
            
            %fba3
            p.formulateFBA([], [], false);
            FbaLPWriter.write_fba3_stoichiometry(sim, [dirName filesep 'metabolism.full.stm'], 'mmol/gDCW/h');
            
            p.formulateFBA([], [], true);
            FbaLPWriter.write_fba3_stoichiometry(sim, [dirName filesep 'metabolism.reduced.stm'], 'mmol/gDCW/h');
            
            FbaLPWriter.write_fba3_deletion(sim, [dirName filesep 'metabolismDeletions.txt']);
        end
        
        %SBML format for use with COBRA toolbox
        %- http://lpsolve.sourceforge.net/5.5/lp-format.htm
        function write_sbml(sim, fileName, units)
            import edu.stanford.covert.cell.sim.util.FbaLPWriter;
            import edu.stanford.covert.util.ConstantUtil;
            
            %parse unit option
            %growth (molecules/cell/s) = growth (mmol/gDCW/h) * biomass scale / flux scale
            if nargin < 3
                units = 'molecules/cell/s';
            elseif ~ismember(units, {'mmol/gDCW/h', 'molecules/cell/s'})
                throw(MException('FbaLPWriter:error', 'Invalid units'))
            end
                        
            [model, compartmentIDs, compartmentNames] = FbaLPWriter.make_cobra_model(sim, units);
            writeCbModel(model, 'sbml', fileName, compartmentIDs, compartmentNames);
            movefile([fileName '.xml'], fileName);
        end

        %MATLAB struct format for use with COBRA toolbox
        %- http://lpsolve.sourceforge.net/5.5/lp-format.htm
        function write_cobra_matlab(sim, fileName, units)
            import edu.stanford.covert.cell.sim.util.FbaLPWriter;
            import edu.stanford.covert.util.ConstantUtil;

            %growth (molecules/cell/s) = growth (mmol/gDCW/h) * biomass scale / flux scale
            if nargin < 3
                units = 'mmol/gDCW/h';
            elseif ~ismember(units, {'mmol/gDCW/h', 'molecules/cell/s'})
                throw(MException('FbaLPWriter:error', 'Invalid units'))
            end

            if ~strcmp(substr(fileName, -4, 4), '.mat')
                fileName = [fileName '.mat'];
            end
            metabolism_full_mmol_gDCW_h = FbaLPWriter.make_cobra_model(sim, units);
            save(fileName, 'metabolism_full_mmol_gDCW_h');
        end
        
        %MATLAB struct format for use with COBRA toolbox
        %- http://lpsolve.sourceforge.net/5.5/lp-format.htm
        function [model, compartmentIDs, compartmentNames] = make_cobra_model(sim, units)
            import edu.stanford.covert.cell.sim.util.FbaLPWriter;
            import edu.stanford.covert.util.ConstantUtil;
            
            %parse unit option
            %growth (molecules/cell/s) = growth (mmol/gDCW/h) * biomass scale / flux scale
            if nargin < 2
                units = 'molecules/cell/s';
            elseif ~ismember(units, {'mmol/gDCW/h', 'molecules/cell/s'})
                throw(MException('FbaLPWriter:error', 'Invalid units'))
            end
            
            %references
            p = sim.process('Metabolism');
            
            %IDs, names
            [compartmentIDs, compartmentNames, substrateIDs, ~, substrateNames, reactionIDs, ~, reactionNames] = ...
                edu.stanford.covert.cell.sim.util.FbaLPWriter.calcIDs(sim, true);
            compartmentIDs = [compartmentIDs; 'none'];
            compartmentNames = [compartmentNames; 'None'];
            
            idx = 0;
            for i = 1:numel(substrateIDs)
                if isempty(substrateIDs{i}) || strcmp(substrateIDs{i}, '[none]')
                    idx = idx + 1;
                    substrateIDs{i} = sprintf('Constraint_%03d[none]', idx);
                    substrateNames{i} = sprintf('Constraint %d', idx);
                end
            end
            
            %scale units
            [fbaReactionStoichiometryMatrix, fluxBounds, fbaObjective] = FbaLPWriter.scaleNetwork(sim, units, ConstantUtil.realmax);
            
            enzGeneComp = p.enzymeGeneComposition();
            tmp = any(enzGeneComp, 2);
            geneIDs = p.gene.wholeCellModelIDs(tmp);
            rxnGeneMat = (p.fbaReactionCatalysisMatrix * enzGeneComp(tmp, :)') ~= 0;            
            grRules = cell(size(reactionIDs));
            for i = 1:numel(reactionIDs)
                grRules{i} = strjoin(' and ', geneIDs{rxnGeneMat(i, :)});
            end
            
            model = struct();
            model.description = 'Mycoplasma genitalium metabolic sub-model';
            model.rxns = reactionIDs;
            model.mets = substrateIDs;
            model.S = fbaReactionStoichiometryMatrix;
            model.lb = fluxBounds(:, 1);
            model.ub = fluxBounds(:, 2);
            model.rev = fluxBounds(:, 1) < 0 & fluxBounds(:, 2) > 0;
            model.c = fbaObjective;
            model.rxnNames = reactionNames;
            model.metNames = substrateNames;
            model.grRules = grRules;
            model.genes = geneIDs;
        end

        %LP format
        %- http://lpsolve.sourceforge.net/5.5/lp-format.htm
        function write_lpsolve_lp(sim, fileName, units)
            import edu.stanford.covert.cell.sim.util.FbaLPWriter;
            import edu.stanford.covert.util.ConstantUtil;
            
            %parse unit option
            %growth (molecules/cell/s) = growth (mmol/gDCW/h) * biomass scale / flux scale
            if nargin < 3
                units = 'molecules/cell/s';
            elseif ~ismember(units, {'mmol/gDCW/h', 'molecules/cell/s'})
                throw(MException('FbaLPWriter:error', 'Invalid units'))
            end
            
            %references
            p = sim.process('Metabolism');
            
            %numbers
            nSubstrates = size(p.fbaReactionStoichiometryMatrix, 1);
            nReactions = size(p.fbaReactionStoichiometryMatrix, 2);
            
            %IDs, names
            [compartmentIDs, compartmentNames, substrateIDs, ~, substrateNames, reactionIDs, ~, reactionNames] = ...
                edu.stanford.covert.cell.sim.util.FbaLPWriter.calcIDs(sim);
            
            %scale units
            [fbaReactionStoichiometryMatrix, fluxBounds, fbaObjective] = FbaLPWriter.scaleNetwork(sim, units, ConstantUtil.realmax);
            
            %initialize
            str = [];
            
            %header
            str = [str sprintf('/*\n')];
            str = [str sprintf('%s\n', 'Mycoplasma genitalium metabolism')];
            str = [str sprintf('Generated %s by edu.stanford.covert.cell.sim.util.FbaLPWriter\n', datestr(now, 31))];
            str = [str sprintf('*/\n\n')];
            
            %objective
            if strcmp(units, 'mmol/gDCW/h')
                str = [str sprintf('// Maximize (mmol/gDCW)\n')];
            else
                str = [str sprintf('// Maximize (molecules/cell)\n')];
            end
            
            idxs = find(fbaObjective);
            tmp = cellfun(@(id, coeff) sprintf('%e %s', coeff, id), reactionIDs(idxs), num2cell(fbaObjective(idxs)), 'UniformOutput', false);
            str = [str sprintf('max: %s;\n', strjoin(' + ', tmp{:}))];
            str = [str sprintf('\n')];
            
            %Constraints
            str = [str sprintf('// Metabolites with compartments\n')];
            for i = 1:numel(compartmentIDs);
                str = [str sprintf('// - %s: %s\n', compartmentIDs{i}, compartmentNames{i})];
            end
            for i = 1:nSubstrates
                idxs = find(fbaReactionStoichiometryMatrix(i, :));
                if isempty(idxs)
                    continue;
                end
                
                tmp = cell(size(idxs));
                for j = 1:numel(idxs)
                    if fbaReactionStoichiometryMatrix(i, idxs(j)) == ceil(fbaReactionStoichiometryMatrix(i, idxs(j)))
                        tmp{j} = sprintf('%d %s', fbaReactionStoichiometryMatrix(i, idxs(j)), reactionIDs{idxs(j)});
                    else
                        tmp{j} = sprintf('%e %s', fbaReactionStoichiometryMatrix(i, idxs(j)), reactionIDs{idxs(j)});
                    end
                end
                str = [str sprintf('%s: %s = 0;\n', substrateIDs{i}, strjoin(' + ', tmp{:}))];
            end
            str = [str sprintf('\n')];
            
            %Bounds
            if strcmp(units, 'mmol/gDCW/h')
                str = [str sprintf('// Bounds (mmol/gDCW/h)\n')];
            else
                str = [str sprintf('// Bounds (1/cell/s)\n')];
            end
            for i = 1:nReactions
                if fluxBounds(i, 1) == ceil(fluxBounds(i, 1))
                    lower = sprintf('%d <= ', fluxBounds(i, 1));
                else
                    lower = sprintf('%e <= ', fluxBounds(i, 1));
                end
                
                if fluxBounds(i, 2) == ceil(fluxBounds(i, 2))
                    upper = sprintf(' <= %d', fluxBounds(i, 2));
                else
                    upper = sprintf(' <= %e', fluxBounds(i, 2));
                end
                
                if fluxBounds(i, 1) == fluxBounds(i, 2)
                    if fluxBounds(i, 1) == ceil(fluxBounds(i, 1))
                        str = [str sprintf('%s = %d;\n', reactionIDs{i}, fluxBounds(i, 1))];
                    else
                        str = [str sprintf('%s = %e;\n', reactionIDs{i}, fluxBounds(i, 1))];
                    end
                else
                    if ~isempty(lower) || ~isempty(upper)
                        str = [str sprintf('%s%s%s;\n', lower, reactionIDs{i}, upper)];
                    end
                end
            end
            str = [str sprintf('\n')];
            
            %define variables
            %str = [str sprintf('// Declare reactions\n')];
            %str = [str sprintf('free %s;\n', strjoin(', ', reactionIDs{:}))];
            
            %write file
            fid = fopen(fileName, 'w');
            fwrite(fid, str);
            fclose(fid);
        end
        
        %LP format
        %- http://lpsolve.sourceforge.net/5.5/CPLEX-format.htm
        %- http://www2.isye.gatech.edu/~wcook/qsopt/hlp/ff_lp_format.htm
        function write_cplex_lp(sim, fileName, units)
            import edu.stanford.covert.cell.sim.util.FbaLPWriter;
            import edu.stanford.covert.util.ConstantUtil;
            
            %parse unit option
            %growth (molecules/cell/s) = growth (mmol/gDCW/h) * biomass scale / flux scale
            if nargin < 3
                units = 'molecules/cell/s';
            elseif ~ismember(units, {'mmol/gDCW/h', 'molecules/cell/s'})
                throw(MException('FbaLPWriter:error', 'Invalid units'))
            end
            
            %references
            p = sim.process('Metabolism');
            
            %numbers
            nSubstrates = size(p.fbaReactionStoichiometryMatrix, 1);
            nReactions = size(p.fbaReactionStoichiometryMatrix, 2);
            
            %IDs, names
            [compartmentIDs, compartmentNames, substrateIDs, ~, ~, reactionIDs, ~, ~] = ...
                edu.stanford.covert.cell.sim.util.FbaLPWriter.calcIDs(sim);
            
            %scale units
            [fbaReactionStoichiometryMatrix, fluxBounds, fbaObjective] = FbaLPWriter.scaleNetwork(sim, units, ConstantUtil.realmax);
            
            %initialize
            str = [];
            
            %header
            str = [str sprintf('\\ %s\n', 'Mycoplasma genitalium metabolism')];
            str = [str sprintf('\\ Generated %s by edu.stanford.covert.cell.sim.util.FbaLPWriter\n\n', datestr(now, 31))];
            
            %objective
            str = [str sprintf('Maximize\n')];
            if strcmp(units, 'mmol/gDCW/h')
                str = [str sprintf('\\ mmol/gDCW\n')];
            else
                str = [str sprintf('\\ molecules/cell\n')];
            end
            
            idxs = find(fbaObjective);
            tmp = cellfun(@(id, coeff) sprintf('%e %s', coeff, id), reactionIDs(idxs), num2cell(fbaObjective(idxs)), 'UniformOutput', false);
            str = [str strrep(sprintf(' obj: %s\n', strjoin(' + ', tmp{:})), ' + -', ' -')];
            
            %Constraints
            str = [str sprintf('Subject To\n')];
            str = [str sprintf('\\ Metabolites\n')];
            str = [str sprintf('\\ Compartments\n')];
            for i = 1:numel(compartmentIDs);
                str = [str sprintf('\\ - %s: %s\n', compartmentIDs{i}, compartmentNames{i})];
            end
            for i = 1:nSubstrates
                idxs = find(fbaReactionStoichiometryMatrix(i, :));
                if isempty(idxs)
                    continue;
                end
                
                tmp = cell(size(idxs));
                for j = 1:numel(idxs)
                    if fbaReactionStoichiometryMatrix(i, idxs(j)) == ceil(fbaReactionStoichiometryMatrix(i, idxs(j)))
                        tmp{j} = sprintf('%d %s', fbaReactionStoichiometryMatrix(i, idxs(j)), reactionIDs{idxs(j)});
                    else
                        tmp{j} = sprintf('%e %s', fbaReactionStoichiometryMatrix(i, idxs(j)), reactionIDs{idxs(j)});
                    end
                end
                str = [str strrep(sprintf(' %s: %s = 0\n', substrateIDs{i}, strjoin(' + ', tmp{:})), ' + -', ' -')];
            end
            
            %Bounds
            str = [str sprintf('Bounds\n')];
            if strcmp(units, 'mmol/gDCW/h')
                str = [str sprintf('\\ mmol/gDCW/h\n')];
            else
                str = [str sprintf('\\ 1/cell/s\n')];
            end
            for i = 1:nReactions
                if fluxBounds(i, 1) == ceil(fluxBounds(i, 1))
                    lower = sprintf('%d', fluxBounds(i, 1));
                else
                    lower = sprintf('%e', fluxBounds(i, 1));
                end
                
                if fluxBounds(i, 2) == ceil(fluxBounds(i, 2))
                    upper = sprintf('%d', fluxBounds(i, 2));
                else
                    upper = sprintf('%e', fluxBounds(i, 2));
                end
                
                if fluxBounds(i, 1) == fluxBounds(i, 2)
                    str = [str sprintf(' %s = %s\n', reactionIDs{i}, lower)];
                else
                    if ~isempty(lower)
                        str = [str sprintf(' %s >= %s\n', reactionIDs{i}, lower)];
                    end
                    if ~isempty(upper)
                        str = [str sprintf(' %s <= %s\n', reactionIDs{i}, upper)];
                    end
                end
            end
            
            %End
            str = [str sprintf('End\n')];
            
            %write file
            fid = fopen(fileName, 'w');
            fwrite(fid, str);
            fclose(fid);
        end
        
        %LP format
        %- http://www.gurobi.com/documentation/5.5/reference-manual/node902
        %- http://www2.isye.gatech.edu/~wcook/qsopt/hlp/ff_lp_format.htm
        function write_gurobi_lp(sim, fileName, units)
            import edu.stanford.covert.cell.sim.util.FbaLPWriter;
            import edu.stanford.covert.util.ConstantUtil;
            
            %parse unit option
            %growth (molecules/cell/s) = growth (mmol/gDCW/h) * biomass scale / flux scale
            if nargin < 3
                units = 'molecules/cell/s';
            elseif ~ismember(units, {'mmol/gDCW/h', 'molecules/cell/s'})
                throw(MException('FbaLPWriter:error', 'Invalid units'))
            end
            
            %references
            p = sim.process('Metabolism');
            
            %numbers
            nSubstrates = size(p.fbaReactionStoichiometryMatrix, 1);
            nReactions = size(p.fbaReactionStoichiometryMatrix, 2);
            
            %IDs, names
            [compartmentIDs, compartmentNames, substrateIDs, ~, ~, reactionIDs, ~, ~] = ...
                edu.stanford.covert.cell.sim.util.FbaLPWriter.calcIDs(sim);
            
            %scale units
            [fbaReactionStoichiometryMatrix, fluxBounds, fbaObjective] = FbaLPWriter.scaleNetwork(sim, units, ConstantUtil.realmax);
            
            %initialize
            str = [];
            
            %header
            str = [str sprintf('\\ %s\n', 'Mycoplasma genitalium metabolism')];
            str = [str sprintf('\\ Generated %s by edu.stanford.covert.cell.sim.util.FbaLPWriter\n\n', datestr(now, 31))];
            
            %objective
            str = [str sprintf('Maximize\n')];
            if strcmp(units, 'mmol/gDCW/h')
                str = [str sprintf('\\ mmol/gDCW\n')];
            else
                str = [str sprintf('\\ molecules/cell\n')];
            end
            
            idxs = find(fbaObjective);
            tmp = cellfun(@(id, coeff) sprintf('%e %s', coeff, id), reactionIDs(idxs), num2cell(fbaObjective(idxs)), 'UniformOutput', false);
            str = [str strrep(sprintf(' obj: %s\n', strjoin(' + ', tmp{:})), ' + -', ' -')];
            
            %Constraints
            str = [str sprintf('Subject To\n')];
            str = [str sprintf('\\ Metabolites\n')];
            str = [str sprintf('\\ Compartments\n')];
            for i = 1:numel(compartmentIDs);
                str = [str sprintf('\\ - %s: %s\n', compartmentIDs{i}, compartmentNames{i})];
            end
            for i = 1:nSubstrates
                idxs = find(fbaReactionStoichiometryMatrix(i, :));
                if isempty(idxs)
                    continue;
                end
                
                tmp = cell(size(idxs));
                for j = 1:numel(idxs)
                    if fbaReactionStoichiometryMatrix(i, idxs(j)) == ceil(fbaReactionStoichiometryMatrix(i, idxs(j)))
                        tmp{j} = sprintf('%d %s', fbaReactionStoichiometryMatrix(i, idxs(j)), reactionIDs{idxs(j)});
                    else
                        tmp{j} = sprintf('%e %s', fbaReactionStoichiometryMatrix(i, idxs(j)), reactionIDs{idxs(j)});
                    end
                end
                str = [str strrep(sprintf(' %s: %s = 0\n', substrateIDs{i}, strjoin(' + ', tmp{:})), ' + -', ' -')];
            end
            
            %Bounds
            str = [str sprintf('Bounds\n')];
            if strcmp(units, 'mmol/gDCW/h')
                str = [str sprintf('\\ mmol/gDCW/h\n')];
            else
                str = [str sprintf('\\ 1/cell/s\n')];
            end
            for i = 1:nReactions
                if fluxBounds(i, 1) == ceil(fluxBounds(i, 1))
                    lower = sprintf('%d', fluxBounds(i, 1));
                else
                    lower = sprintf('%e', fluxBounds(i, 1));
                end
                
                if fluxBounds(i, 2) == ceil(fluxBounds(i, 2))
                    upper = sprintf('%d', fluxBounds(i, 2));
                else
                    upper = sprintf('%e', fluxBounds(i, 2));
                end
                
                if fluxBounds(i, 1) == fluxBounds(i, 2)
                    str = [str sprintf('%s = %s\n', reactionIDs{i}, upper)];
                else
                    str = [str sprintf('%s <= %s <= %s\n', lower, reactionIDs{i}, upper)];
                end
            end
            
            %End
            str = [str sprintf('End\n')];
            
            %write file
            fid = fopen(fileName, 'w');
            fwrite(fid, str);
            fclose(fid);
        end
        
        function write_lindo(sim, fileName, units)
            import edu.stanford.covert.cell.sim.util.FbaLPWriter;
            import edu.stanford.covert.util.ConstantUtil;
            
            %parse unit option
            %growth (molecules/cell/s) = growth (mmol/gDCW/h) * biomass scale / flux scale
            if nargin < 3
                units = 'molecules/cell/s';
            elseif ~ismember(units, {'mmol/gDCW/h', 'molecules/cell/s'})
                throw(MException('FbaLPWriter:error', 'Invalid units'))
            end
            
            %references
            p = sim.process('Metabolism');
            
            %numbers
            nSubstrates = size(p.fbaReactionStoichiometryMatrix, 1);
            nReactions = size(p.fbaReactionStoichiometryMatrix, 2);
            
            %IDs, names
            [compartmentIDs, compartmentNames, substrateIDs, ~, substrateNames, reactionIDs, ~, reactionNames] = ...
                edu.stanford.covert.cell.sim.util.FbaLPWriter.calcIDs(sim);
            
            %scale units
            [fbaReactionStoichiometryMatrix, fluxBounds, fbaObjective] = FbaLPWriter.scaleNetwork(sim, units, ConstantUtil.realmax);
            
            %initialize
            str = [];
            
            %header
            str = [str sprintf('! Autogenerated %s by edu.stanford.covert.cell.sim.util.FbaLPWriter\n', datestr(now, 31))];
            str = [str sprintf('TITLE %s\n\n', 'Metabolism')];
            
            if strcmp(units, 'mmol/gDCW/h')
                str = [str sprintf('! Objective function (mmol/gDCW)\n')];
            else
                str = [str sprintf('! Objective function (molecules/cell)\n')];
            end
            str = [str sprintf('MAXIMIZE\n')];
            idxs = find(fbaObjective);
            for i = 1:numel(idxs)
                str = [str sprintf('\t%e %s\n', fbaObjective(idxs(i)), reactionIDs{idxs(i)})];
            end
            str = [str sprintf('\n')];
            
            str = [str sprintf('! Reactions\nSUBJECT TO\n')];
            
            for i = 1:nSubstrates
                rxnIdxs = find(p.fbaReactionStoichiometryMatrix(i, :));
                if isempty(rxnIdxs)
                    continue;
                end
                
                constraint = [];
                for j = 1:numel(rxnIdxs)
                    coeff = fbaReactionStoichiometryMatrix(i, rxnIdxs(j));
                    if coeff == floor(coeff)
                        constraint = [constraint sprintf('%d %s ', coeff, reactionIDs{rxnIdxs(j)})];
                    else
                        constraint = [constraint sprintf('%e %s ', coeff, reactionIDs{rxnIdxs(j)})];
                    end
                end
                str = [str sprintf(' %s) %s= 0 ! %s\n', substrateIDs{i}, constraint, substrateNames{i})];
            end
            
            str = [str sprintf('\nEND\n\n')];
            
            if strcmp(units, 'mmol/gDCW/h')
                str = [str sprintf('! Reaction Bounds (mmol/gDCW/h)\n')];
            else
                str = [str sprintf('! Reaction Bounds (1/cell/s)\n')];
            end
            
            for i = 1:nReactions
                enzIdx = find(p.fbaReactionCatalysisMatrix(i, :));
                
                if fluxBounds(i, 1) == floor(fluxBounds(i, 1))
                    lb = sprintf('%d', fluxBounds(i, 1));
                else
                    lb = sprintf('%d', fluxBounds(i, 1));
                end
                if fluxBounds(i, 2) == floor(fluxBounds(i, 2))
                    ub = sprintf('%d', fluxBounds(i, 2));
                else
                    ub = sprintf('%d', fluxBounds(i, 2));
                end
                
                if ~isempty(enzIdx)
                    str = [str sprintf('SLB %s %s ! %s (%s)\n', reactionIDs{i}, lb, reactionNames{i}, p.enzymeWholeCellModelIDs{enzIdx})];
                    str = [str sprintf('SUB %s %s ! %s (%s)\n', reactionIDs{i}, ub, reactionNames{i}, p.enzymeWholeCellModelIDs{enzIdx})];
                else
                    str = [str sprintf('SLB %s %s ! %s\n', reactionIDs{i}, lb, reactionNames{i})];
                    str = [str sprintf('SUB %s %s ! %s\n', reactionIDs{i}, ub, reactionNames{i})];
                end
            end
            
            str = [str sprintf('\n')];
            str = [str sprintf('! Compartments\n')];
            for i = 1:numel(compartmentIDs);
                str = [str sprintf('! %s\t%s\n', compartmentIDs{i}, compartmentNames{i})];
            end
            
            fid = fopen(fileName, 'w');
            fwrite(fid, str);
            fclose(fid);
        end
        
        function write_fba3_stoichiometry(sim, fileName, units)
            import edu.stanford.covert.cell.sim.util.FbaLPWriter;
            import edu.stanford.covert.util.ConstantUtil;
            
            %references
            p = sim.process('Metabolism');
            
            %numbers
            nSubstrates = size(p.reactionStoichiometryMatrix, 1);
            nCompartments = size(p.reactionStoichiometryMatrix, 3);
            
            %ids
            [compartmentIDs, compartmentNames, ~, substrateIDs, substrateNames, ~, reactionIDs, reactionNames, reversible] = ...
                edu.stanford.covert.cell.sim.util.FbaLPWriter.calcIDs(sim);
            
            %scale units
            [fbaSMat, fluxBounds, fbaObj] = FbaLPWriter.scaleNetwork(sim, units, 1000);
            fbaObj = fbaObj * 10^max(0, -floor(log10(min(abs(fbaObj(fbaObj ~= 0))))));
            
            %initialize
            str = [];
            
            %header
            str = [str sprintf('// %s //\n', repmat('*', 1, 120))];
            str = [str sprintf('// Metabolic network //\n')];
            str = [str sprintf('// Autogenerated %s by edu.stanford.covert.cell.sim.util.FbaLPWriter //\n', datestr(now, 31))];
            str = [str sprintf('// %s //\n', repmat('*', 1, 120))];
            str = [str sprintf('\n')];
            str = [str sprintf('\n')];
            
            %compartments
            idWidth = max(cellfun(@length, compartmentIDs));
            nameWidth = max(cellfun(@length, compartmentNames));
            
            str = [str sprintf('// %s //\n', repmat('*', 1, 120))];
            str = [str sprintf('// Compartments %s //\n', repmat(' ', 1, idWidth + nameWidth + 3 - numel('Compartments')))];
            for i = 1:nCompartments
                str = [str sprintf('// %s\t%s //\n', ...
                    [compartmentIDs{i} repmat(' ', 1, idWidth - numel(compartmentIDs{i}))], ...
                    [compartmentNames{i} repmat(' ', 1, nameWidth - numel(compartmentNames{i}))])];
            end
            str = [str sprintf('// %s //\n', repmat('*', 1, 120))];
            str = [str sprintf('\n')];
            str = [str sprintf('\n')];
            
            %substrates
            idWidth = max(cellfun(@length, substrateIDs));
            nameWidth = max(cellfun(@length, substrateNames));
            
            [subIdxs, ~] = ind2sub([nSubstrates, nCompartments], p.substrateIndexs_fba);
            
            str = [str sprintf('// %s //\n', repmat('*', 1, 120))];
            str = [str sprintf('// Substrates %s //\n', repmat(' ', 1, idWidth + nameWidth +5 - length('Substrates')))];
            str = [str sprintf('// %s //\n', repmat(' ', 1, idWidth + nameWidth + 6))];
            
            str = [str sprintf('// \tMetabolites %s //\n', repmat(' ', 1, idWidth + nameWidth +4 - length('Metabolites')))];
            for i = 1:size(p.substrateMetaboliteLocalIndexs, 1)
                idx = p.fbaSubstrateIndexs_substrates(p.substrateMetaboliteLocalIndexs(i) == subIdxs);
                if isempty(idx)
                    continue;
                elseif numel(idx)==1
                    id = substrateIDs{idx};
                else
                    idx = idx(1);
                    id = substrateIDs{idx}(1:end-1);
                end
                name = substrateNames{idx};
                str = [str sprintf('// \t\t%s\t%s //\n', ...
                    [id repmat(' ', 1, idWidth - numel(id))], ...
                    [name repmat(' ', 1, nameWidth - numel(name))])]; %#ok<*AGROW>
            end
            str = [str sprintf('// %s //\n', repmat(' ', 1, idWidth + nameWidth + 6))];
            
            str = [str sprintf('// \tProtein Monomers %s //\n', repmat(' ', 1, idWidth + nameWidth +4 - length('Protein Monomers')))];
            for i = 1:size(p.substrateMonomerLocalIndexs, 1)
                idx = p.fbaSubstrateIndexs_substrates(p.substrateMonomerLocalIndexs(i) == subIdxs);
                if isempty(idx)
                    continue;
                elseif numel(idx)==1
                    id = substrateIDs{idx};
                else
                    idx = idx(1);
                    id = substrateIDs{idx}(1:end-1);
                end
                name = substrateNames{idx};
                str = [str sprintf('// \t\t%s\t%s //\n', ...
                    [id repmat(' ', 1, idWidth - numel(id))], ...
                    [name repmat(' ', 1, nameWidth - numel(name))])]; %#ok<*AGROW>
            end
            str = [str sprintf('// %s //\n', repmat(' ', 1, idWidth + nameWidth + 6))];
            
            str = [str sprintf('// \tProtein Complexs %s //\n', repmat(' ', 1, idWidth + nameWidth +4 - length('Protein Complexs')))];
            for i = 1:size(p.substrateComplexLocalIndexs, 1)
                idx = p.fbaSubstrateIndexs_substrates(p.substrateComplexLocalIndexs(i) == subIdxs);
                if isempty(idx)
                    continue;
                elseif numel(idx)==1
                    id = substrateIDs{idx};
                else
                    idx = idx(1);
                    id = substrateIDs{idx}(1:end-1);
                end
                name = substrateNames{idx};
                str = [str sprintf('// \t\t%s\t%s //\n', ...
                    [id repmat(' ', 1, idWidth - numel(id))], ...
                    [name repmat(' ', 1, nameWidth - numel(name))])]; %#ok<*AGROW>
            end
            str = [str sprintf('// %s //\n', repmat(' ', 1, idWidth + nameWidth + 6))];
            
            str = [str sprintf('// \tInternal exchange constraints %s //\n', repmat(' ', 1, idWidth + nameWidth +4 - length('Internal exchange constraints')))];
            for i = 1:size(p.fbaSubstrateIndexs_metaboliteInternalExchangeConstraints, 1)
                idx = p.fbaSubstrateIndexs_metaboliteInternalExchangeConstraints(i);
                if isempty(idx)
                    continue;
                elseif numel(idx)==1
                    id = substrateIDs{idx};
                else
                    idx = idx(1);
                    id = substrateIDs{idx}(1:end-1);
                end
                name = substrateNames{idx};
                str = [str sprintf('// \t\t%s\t%s //\n', ...
                    [id repmat(' ', 1, idWidth - numel(id))], ...
                    [name repmat(' ', 1, nameWidth - numel(name))])]; %#ok<*AGROW>
            end
            str = [str sprintf('// %s //\n', repmat(' ', 1, idWidth + nameWidth + 6))];
            
            str = [str sprintf('// \tPseudospecies %s //\n', repmat(' ', 1, idWidth + nameWidth +4 - length('Pseudospecies')))];
            str = [str sprintf('// \t\t%s\t%s //\n', ...
                ['Biomass' repmat(' ', 1, idWidth - numel('Biomass'))], ...
                ['Biomass' repmat(' ', 1, nameWidth - numel('Biomass'))])]; %#ok<*AGROW>
            str = [str sprintf('// %s //\n', repmat('*', 1, 120))];
            str = [str sprintf('\n')];
            str = [str sprintf('\n')];
            
            %reactions
            str = [str sprintf('// %s //\n', repmat('*', 1, 120))];
            str = [str sprintf('// Reactions //\n')];
            str = [str sprintf('\n')];
            
            %biomass production in mmol/gDCW
            str = [str sprintf('// \tBiomass production //\n')];
            
            biomassScaleFactor = 100;
            idx = p.fbaReactionIndexs_biomassProduction;
            stoichiometry = [];
            stoichiometry2 = [];
            
            scIdxs = p.fbaSubstrateIndexs_substrates(fbaSMat(p.fbaSubstrateIndexs_substrates, idx) < 0);
            for j = 1:numel(scIdxs)
                val = fbaSMat(scIdxs(j), idx);
                
                stoichiometry = [stoichiometry sprintf('%s%f %s ', repmat('+', 1, val >= 0), val * biomassScaleFactor, substrateIDs{scIdxs(j)})];
                stoichiometry2 = [stoichiometry2 sprintf('// %17.6f %s%s //\n', ...
                    val, ...
                    substrateIDs{scIdxs(j)}, ...
                    repmat(' ', 1, idWidth-numel(substrateIDs{scIdxs(j)})))];
            end
            
            scIdxs = p.fbaSubstrateIndexs_substrates(fbaSMat(p.fbaSubstrateIndexs_substrates, idx) > 0);
            for j = 1:numel(scIdxs)
                val = fbaSMat(scIdxs(j), idx);
                
                stoichiometry = [stoichiometry sprintf('%s%f %s ', repmat('+', 1, val >= 0), val * biomassScaleFactor, substrateIDs{scIdxs(j)})];
                stoichiometry2 = [stoichiometry2 sprintf('// %17.6f %s%s //\n', ...
                    val, ...
                    substrateIDs{scIdxs(j)}, ...
                    repmat(' ', 1, idWidth-numel(substrateIDs{scIdxs(j)})))];
            end
            stoichiometry = stoichiometry(1:end-1);
            
            str = [str stoichiometry2];
            str = [str sprintf('\n')];
            
            str = [str sprintf('// Biomass production in mmol/gDCW //\n')];
            spaces = find(stoichiometry == ' ');
            spaces = [spaces(2:2:end) length(stoichiometry)+1];
            nSplits = ceil((length(spaces)-1)/10);
            splits = [0 spaces(ceil((1:nSplits)/nSplits * length(spaces)))];
            assertEqual(splits, unique(splits));
            growthStoichiometry = [];
            for j = 1:nSplits
                str = [str sprintf('%s +%d SbBmss%02d 0 SbGrth%02d // Sub growth - %d //\n', stoichiometry(splits(j)+1:splits(j+1)-1), biomassScaleFactor, j, j, j)];
                growthStoichiometry = [growthStoichiometry sprintf('-1 SbBmss%02d ', j)];
            end
            
            str = [str sprintf('%s +1 Biomass 0 %s // %s //\n', growthStoichiometry(1:end-1), reactionIDs{idx}, reactionNames{idx})];
            str = [str sprintf('\n')];
            
            str = [str sprintf('// \tMetabolic Reactions //\n')];
            for i = 1:numel(p.reactionIndexs_chemical)
                idx = p.fbaReactionIndexs_metabolicConversion(p.reactionIndexs_chemical(i) == p.reactionIndexs_fba);
                if isempty(idx)
                    continue;
                end
                stoichiometry = [];
                
                scIdxs = find(fbaSMat(:, idx) < 0);
                for j = 1:numel(scIdxs)
                    val = fbaSMat(scIdxs(j), idx);
                    stoichiometry = [stoichiometry sprintf('%s%d %s ', repmat('+', 1, val>=0), val, substrateIDs{scIdxs(j)})];
                end
                
                scIdxs = find(fbaSMat(:, idx) > 0);
                for j = 1:numel(scIdxs)
                    val = fbaSMat(scIdxs(j), idx);
                    stoichiometry = [stoichiometry sprintf('%s%d %s ', repmat('+', 1, val>=0), val, substrateIDs{scIdxs(j)})];
                end
                
                eIdx = find(p.fbaReactionCatalysisMatrix(idx, :));
                
                if ~isempty(eIdx)
                    enzymeID = p.enzymeWholeCellModelIDs{eIdx};
                    enzymeName = p.enzymeNames{eIdx};
                    str = [str sprintf('%s0 %s // %s (%s, %s) //\n', stoichiometry, reactionIDs{idx}, reactionNames{idx}, enzymeID, enzymeName)];
                else
                    str = [str sprintf('%s0 %s // %s //\n', stoichiometry, reactionIDs{idx}, reactionNames{idx})];
                end
            end
            str = [str sprintf('\n')];
            
            str = [str sprintf('// \tTransport Reactions //\n')];
            for i = 1:numel(p.reactionIndexs_transport)
                idx = p.fbaReactionIndexs_metabolicConversion(p.reactionIndexs_transport(i) == p.reactionIndexs_fba);
                if isempty(idx)
                    continue;
                end
                stoichiometry = [];
                
                scIdxs = find(fbaSMat(:, idx) < 0);
                for j = 1:numel(scIdxs)
                    val = fbaSMat(scIdxs(j), idx);
                    stoichiometry = [stoichiometry sprintf('%s%d %s ', repmat('+', 1, val>=0), val, substrateIDs{scIdxs(j)})];
                end
                
                scIdxs = find(fbaSMat(:, idx) > 0);
                for j = 1:numel(scIdxs)
                    val = fbaSMat(scIdxs(j), idx);
                    stoichiometry = [stoichiometry sprintf('%s%d %s ', repmat('+', 1, val>=0), val, substrateIDs{scIdxs(j)})];
                end
                
                eIdx = find(p.fbaReactionCatalysisMatrix(idx, :));
                if ~isempty(eIdx)
                    enzymeID = p.enzymeWholeCellModelIDs{eIdx};
                    enzymeName = p.enzymeNames{eIdx};
                    str = [str sprintf('%s0 %s // %s (%s, %s) //\n', stoichiometry, reactionIDs{idx}, reactionNames{idx}, enzymeID, enzymeName)];
                else
                    str = [str sprintf('%s0 %s // %s //\n', stoichiometry, reactionIDs{idx}, reactionNames{idx})];
                end
            end
            str = [str sprintf('\n')];
            
            str = [str sprintf('// \tMetabolite External Exchange Reactions //\n')];
            for i = 1:numel(p.fbaReactionIndexs_metaboliteExternalExchange)
                idx = p.fbaReactionIndexs_metaboliteExternalExchange(i);
                scIdxs = find(fbaSMat(:, idx));
                
                stoichiometry = [];
                for j = 1:numel(scIdxs)
                    val = fbaSMat(scIdxs(j), idx);
                    stoichiometry = [stoichiometry sprintf('%s%d %s ', repmat('+', 1, val>=0), val, substrateIDs{scIdxs(j)})];
                end
                str = [str sprintf('%s0 %s // %s //\n', stoichiometry, reactionIDs{idx}, reactionNames{idx})];
            end
            str = [str sprintf('\n')];
            
            str = [str sprintf('// \tMetabolite Internal Exchange Reactions //\n')];
            for i = 1:numel(p.fbaReactionIndexs_metaboliteInternalExchange)
                idx = p.fbaReactionIndexs_metaboliteInternalExchange(i);
                if isempty(idx)
                    continue;
                end
                stoichiometry = [];
                
                scIdxs = find(fbaSMat(:, idx) < 0);
                for j = 1:numel(scIdxs)
                    val = fbaSMat(scIdxs(j), idx);
                    stoichiometry = [stoichiometry sprintf('%s%d %s ', repmat('+', 1, val>=0), val, substrateIDs{scIdxs(j)})];
                end
                
                scIdxs = find(fbaSMat(:, idx) > 0);
                for j = 1:numel(scIdxs)
                    val = fbaSMat(scIdxs(j), idx);
                    stoichiometry = [stoichiometry sprintf('%s%d %s ', repmat('+', 1, val>=0), val, substrateIDs{scIdxs(j)})];
                end
                
                eIdx = find(p.fbaReactionCatalysisMatrix(idx, :));
                if ~isempty(eIdx)
                    enzymeID = p.enzymeWholeCellModelIDs{eIdx};
                    enzymeName = p.enzymeNames{eIdx};
                    str = [str sprintf('%s0 %s // %s (%s, %s) //\n', stoichiometry, reactionIDs{idx}, reactionNames{idx}, enzymeID, enzymeName)];
                else
                    str = [str sprintf('%s0 %s // %s //\n', stoichiometry, reactionIDs{idx}, reactionNames{idx})];
                end
            end
            str = [str sprintf('\n')];
            
            str = [str sprintf('// \tBiomass Exchange Reactions //\n')];
            idx = p.fbaReactionIndexs_biomassExchange;
            scIdxs = find(fbaSMat(:, idx));
            stoichiometry = [];
            for j = 1:numel(scIdxs)
                val = fbaSMat(scIdxs(j), idx);
                stoichiometry = [stoichiometry sprintf('%s%d %s ', repmat('+', 1, val>=0), val, substrateIDs{scIdxs(j)})];
            end
            str = [str sprintf('%s0 %s // %s //\n', stoichiometry, reactionIDs{idx}, reactionNames{idx})];
            str = [str sprintf('// %s //\n', repmat('*', 1, 120))];
            str = [str sprintf('\n')];
            str = [str sprintf('\n')];
            
            %objective
            str = [str sprintf('// %s //\n', repmat('*', 1, 120))];
            str = [str sprintf('// Objective //\n')];
            str = [str sprintf('0.0 end\n')];
            str = [str sprintf('\n')];
            str = [str sprintf('end E 0\n')];
            str = [str sprintf('max\n')];
            idxs = find(fbaObj);
            for i = 1:numel(idxs)
                str = [str sprintf('%.6f %s\n', fbaObj(idxs(i)), reactionIDs{idxs(i)})];
            end
            str = [str sprintf('0 end\n')];
            str = [str sprintf('// %s //\n', repmat('*', 1, 120))];
            str = [str sprintf('\n')];
            str = [str sprintf('\n')];
            
            %reaction bounds
            str = [str sprintf('// %s //\n', repmat('*', 1, 120))];
            str = [str sprintf('// Reaction Bounds (mmol/gDCW/h) //\n')];
            for i = 1:size(fluxBounds, 1)
                if ceil(fluxBounds(i, 1)) == fluxBounds(i, 1)
                    loStr = sprintf('%d', fluxBounds(i, 1));
                else
                    loStr = sprintf('%.6f', fluxBounds(i, 1));
                end
                if ceil(fluxBounds(i, 2)) == fluxBounds(i, 2)
                    hiStr = sprintf('%d', fluxBounds(i, 2));
                else
                    hiStr = sprintf('%.6f', fluxBounds(i, 2));
                end
                str = [str sprintf('%s %s %s\n', loStr, reactionIDs{i}, hiStr)];
                % if reversible(i) && (fbaReactionBounds(i, 1) ~= -1000 || fbaReactionBounds(i, 2) ~= 1000)
                %    str = [str sprintf('%.6f %s %.6f\n', fbaReactionBounds(i, 1), reactionIDs{i}, fbaReactionBounds(i, 2))];
                % elseif ~reversible(i) && (fbaReactionBounds(i, 1) ~= 0 || fbaReactionBounds(i, 2) ~= 1000)
                %    str = [str sprintf('%.6f %s %.6f\n', fbaReactionBounds(i, 1), reactionIDs{i}, fbaReactionBounds(i, 2))];
                % end
            end
            str = [str sprintf('\n')];
            
            str = [str sprintf('// End of reaction Bounds //\n')];
            str = [str sprintf('0 end 0\n')];
            str = [str sprintf('// %s //\n', repmat('*', 1, 120))];
            
            fid = fopen(fileName, 'w');
            fwrite(fid, str);
            fclose(fid);
        end
        
        function write_fba3_deletion(sim, fileName)
            [~, ~, ~, ~, ~, ~, reactionIDs] = edu.stanford.covert.cell.sim.util.FbaLPWriter.calcIDs(sim);
            for i = 1:numel(reactionIDs)
                reactionIDs{i} = regexprep(reactionIDs{i}, 'R$', '');
                reactionIDs{i} = regexprep(reactionIDs{i}, '\d+$', '');
            end
            
            %initialize file
            str = [];
            
            %genes - reactions
            m = sim.process('Metabolism');
            geneEnzComp = m.enzymeGeneComposition ~= 0;
            geneSubComp = m.substrateGeneComposition ~= 0;
            
            [sIdxs, ~] = ind2sub([size(m.reactionStoichiometryMatrix, 1) size(m.reactionStoichiometryMatrix, 3)], m.substrateIndexs_fba);
            for i = 1:numel(sim.gene.wholeCellModelIDs)
                enzIdxs = find(geneEnzComp(i, :));
                subIdxs = find(geneSubComp(i, :));
                
                if isempty(enzIdxs) && isempty(subIdxs)
                    continue;
                end
                
                fbaSubIdxs = m.fbaSubstrateIndexs_substrates(ismember(sIdxs, subIdxs));
                rxnIdxs = any(m.fbaReactionCatalysisMatrix(:, enzIdxs), 2) | any(m.fbaReactionStoichiometryMatrix(fbaSubIdxs, :), 1)';
                rxnIDs = unique(reactionIDs(rxnIdxs));
                for j = 1:7:numel(rxnIDs)
                    tmpRxnIDs = rxnIDs(j:min(j+6, end));
                    str = [str sprintf('%d\t%sg\t%s\n', numel(tmpRxnIDs)+1, strrep(sim.gene.wholeCellModelIDs{i},'_',''), strjoin(' ', tmpRxnIDs{:}))];
                end
            end
            
            %terminator
            str = [str sprintf('\n999 end\n')];
            
            %write file
            fid = fopen(fileName, 'w');
            fwrite(fid, str);
            fclose(fid);
        end
        
        function [sMat, bounds, obj] = scaleNetwork(sim, units, maxFlux)
            import edu.stanford.covert.util.ConstantUtil;
            
            p = sim.process('Metabolism');
            mass = sim.state('Mass');
            
            sMat = p.fbaReactionStoichiometryMatrix;
            bounds = p.calcFluxBounds(p.substrates, p.enzymes, p.fbaReactionBounds, p.fbaEnzymeBounds);
            obj = p.fbaObjective;
            if strcmp(units, 'mmol/gDCW/h')
                %biomass composition
                biomassScaleFactor = 1 / (1e-3 * ConstantUtil.nAvogadro) / mass.cellInitialDryWeight;
                sMat(:, p.fbaReactionIndexs_biomassProduction) = ...
                    sMat(:, p.fbaReactionIndexs_biomassProduction) * biomassScaleFactor;
                sMat(p.fbaSubstrateIndexs_biomass, p.fbaReactionIndexs_biomassProduction) = 1;
                
                %fluxes
                fluxScaleFactor = 1 / (ConstantUtil.nAvogadro/1000) * ConstantUtil.secondsPerHour / ...
                    sum(mass.cellDry);
                bounds = bounds * fluxScaleFactor;
                
                %objective
                obj = p.fbaObjective / fluxScaleFactor;
                obj(p.fbaReactionIndexs_biomassProduction) = ...
                    p.fbaObjective(p.fbaReactionIndexs_biomassProduction) / ConstantUtil.secondsPerHour;
            end
            
            bounds(:, 1) = max(-maxFlux, bounds(:, 1));
            bounds(:, 2) = min( maxFlux, bounds(:, 2));
        end
        
        function [compartmentIDs, compartmentNames, ...
                substrateLongIDs, substrateIDs, substrateNames, ...
                reactionLongIDs, reactionIDs, reactionNames, reversible] = ...
                calcIDs(sim, formatSubstrateIdsForSbml)
            if nargin < 2
                formatSubstrateIdsForSbml = false;
            end
                
            p = sim.process('Metabolism');
            
            %% compartments
            compartmentIDs = p.compartment.wholeCellModelIDs(p.substrateMetaboliteCompartmentIndexs(1, :));
            compartmentNames = p.compartment.names(p.substrateMetaboliteCompartmentIndexs(1, :));
            
            %% substrates
            tmp_substrateShortLongIDs = struct;
            tmp_substrateLongIDs = p.substrateWholeCellModelIDs;
            tmp_substrateIDs = cellfun(@(x) strrep(x, '_', ''), p.substrateWholeCellModelIDs, 'UniformOutput', false);
            tmp_substrateNames = p.substrateNames;
            for i = 1:numel(tmp_substrateIDs)
                if numel(tmp_substrateIDs{i}) > 6
                    if ~isfield(tmp_substrateShortLongIDs, upper(tmp_substrateIDs{i}(1:4)))
                        tmp_substrateShortLongIDs.(upper(tmp_substrateIDs{i}(1:4))) = 0;
                    end
                    tmp_substrateShortLongIDs.(upper(tmp_substrateIDs{i}(1:4))) = ...
                        tmp_substrateShortLongIDs.(upper(tmp_substrateIDs{i}(1:4))) + 1;
                    tmp_substrateIDs{i} = sprintf('%s%02d', tmp_substrateIDs{i}(1:4), tmp_substrateShortLongIDs.(upper(tmp_substrateIDs{i}(1:4))));
                end
            end
            assertEqual(numel(tmp_substrateIDs), numel(unique(tmp_substrateIDs)));
            assertIn(max(cellfun(@length, tmp_substrateIDs)), [0 6]);
            
            substrateLongIDs = cell(size(p.fbaReactionStoichiometryMatrix, 1), 1);
            substrateIDs = cell(size(p.fbaReactionStoichiometryMatrix, 1), 1);
            substrateNames = cell(size(p.fbaReactionStoichiometryMatrix, 1), 1);
            substrateCompartmentIDs = cell(size(p.fbaReactionStoichiometryMatrix, 1), 1);
            
            [subIdxs, cmpIdxs] = ind2sub([numel(tmp_substrateIDs) numel(compartmentIDs)], p.substrateIndexs_fba);
            substrateLongIDs(p.fbaSubstrateIndexs_substrates) = tmp_substrateLongIDs(subIdxs);
            substrateIDs(p.fbaSubstrateIndexs_substrates) = tmp_substrateIDs(subIdxs);
            substrateNames(p.fbaSubstrateIndexs_substrates) = tmp_substrateNames(subIdxs);
            cnts = histc(subIdxs, 1:numel(tmp_substrateIDs));
            substrateCompartmentIDs(p.fbaSubstrateIndexs_substrates(cnts(subIdxs) > 1)) = compartmentIDs(cmpIdxs(cnts(subIdxs) > 1));
            
            substrateLongIDs(p.fbaSubstrateIndexs_metaboliteInternalExchangeConstraints) = ...
                cellfun(@(idx) sprintf('Internal_Exchange_Constraint_%02d', idx), ...
                num2cell(1:numel(p.fbaSubstrateIndexs_metaboliteInternalExchangeConstraints))', 'UniformOutput', false);
            substrateIDs(p.fbaSubstrateIndexs_metaboliteInternalExchangeConstraints) = ...
                cellfun(@(idx) sprintf('ExCt%02d', idx), ...
                num2cell(1:numel(p.fbaSubstrateIndexs_metaboliteInternalExchangeConstraints))', 'UniformOutput', false);
            substrateNames(p.fbaSubstrateIndexs_metaboliteInternalExchangeConstraints) = ...
                cellfun(@(idx) ['Internal Exchange Constraint - #' num2str(idx)], ...
                num2cell(1:numel(p.fbaSubstrateIndexs_metaboliteInternalExchangeConstraints))', 'UniformOutput', false);
            substrateCompartmentIDs(p.fbaSubstrateIndexs_metaboliteInternalExchangeConstraints) = {''};
            
            substrateLongIDs{p.fbaSubstrateIndexs_biomass} = 'Biomass';
            substrateIDs{p.fbaSubstrateIndexs_biomass} = 'Biomass';
            substrateNames{p.fbaSubstrateIndexs_biomass} = 'Biomass';
            substrateCompartmentIDs{p.fbaSubstrateIndexs_biomass} = '';
            
            assertIn(max(cellfun(@length, substrateIDs)), [0 8]);
            
            %% reactions
            tmp_reversible = false(size(p.reactionStoichiometryMatrix, 2), 1);
            for i = 1:size(p.reactionStoichiometryMatrix, 2)
                if all(p.reactionBounds(i, :))
                    tmp_reversible(i) = true;
                end
            end
            
            tmp_reactionShortLongIDs = struct;
            tmp_reactionLongIDs = p.reactionWholeCellModelIDs;
            tmp_reactionIDs = cellfun(@(x) strrep(x, '_', ''), tmp_reactionLongIDs, 'UniformOutput', false);
            tmp_reactionNames = p.reactionNames;
            for i = 1:numel(tmp_reactionIDs)
                if ~tmp_reversible(i) && numel(tmp_reactionIDs{i}) > 8
                    if ~isfield(tmp_reactionShortLongIDs, upper(tmp_reactionIDs{i}(1:6)))
                        tmp_reactionShortLongIDs.(upper(tmp_reactionIDs{i}(1:6))) = 0;
                    end
                    tmp_reactionShortLongIDs.(upper(tmp_reactionIDs{i}(1:6))) = ...
                        tmp_reactionShortLongIDs.(upper(tmp_reactionIDs{i}(1:6))) + 1;
                    tmp_reactionIDs{i} = sprintf('%s%02d', tmp_reactionIDs{i}(1:6), tmp_reactionShortLongIDs.(upper(tmp_reactionIDs{i}(1:6))));
                elseif tmp_reversible(i) && numel(tmp_reactionIDs{i}) > 7
                    if ~isfield(tmp_reactionShortLongIDs, upper(tmp_reactionIDs{i}(1:5)))
                        tmp_reactionShortLongIDs.(upper(tmp_reactionIDs{i}(1:5))) = 0;
                    end
                    tmp_reactionShortLongIDs.(upper(tmp_reactionIDs{i}(1:5))) = ...
                        tmp_reactionShortLongIDs.(upper(tmp_reactionIDs{i}(1:5))) + 1;
                    tmp_reactionIDs{i} = sprintf('%s%02dR', tmp_reactionIDs{i}(1:5), tmp_reactionShortLongIDs.(upper(tmp_reactionIDs{i}(1:5))));
                elseif tmp_reversible(i)
                    tmp_reactionIDs{i} = sprintf('%sR', tmp_reactionIDs{i});
                elseif ~tmp_reversible(i) && upper(tmp_reactionIDs{i}(end)) == 'R'
                    if ~isfield(tmp_reactionShortLongIDs, upper(tmp_reactionIDs{i}(1:min(7,end))))
                        tmp_reactionShortLongIDs.(upper(tmp_reactionIDs{i}(1:min(7,end)))) = 0;
                    end
                    tmp_reactionShortLongIDs.(upper(tmp_reactionIDs{i}(1:min(7,end)))) = ...
                        tmp_reactionShortLongIDs.(upper(tmp_reactionIDs{i}(1:min(7,end)))) + 1;
                    tmp_reactionIDs{i} = sprintf('%sF', tmp_reactionIDs{i}(1:min(7,end)));
                end
            end
            
            externalExchangeReactionLongIDs = {};
            externalExchangeReactionIDs = {};
            externalExchangeReactionNames = {};
            for i = 1:numel(p.fbaReactionIndexs_metaboliteExternalExchange)
                sIdx = find(p.fbaReactionStoichiometryMatrix(:, p.fbaReactionIndexs_metaboliteExternalExchange(i)) ~= 0);
                longid = substrateLongIDs{sIdx};
                id = substrateIDs{sIdx};
                name = substrateNames{sIdx};
                externalExchangeReactionLongIDs = [externalExchangeReactionLongIDs; {[longid 'ex']}];
                externalExchangeReactionIDs = [externalExchangeReactionIDs; {[id 'ex']}];
                externalExchangeReactionNames = [externalExchangeReactionNames; {[name ' external exchange']}];
            end
            
            internalExchangeReactionLongIDs = {};
            internalExchangeReactionIDs = {};
            internalExchangeReactionNames = {};
            for i = 1:numel(p.fbaReactionIndexs_metaboliteInternalExchange)
                sIdx = find(p.fbaReactionStoichiometryMatrix(:, p.fbaReactionIndexs_metaboliteInternalExchange(i)) == -1);
                longid = substrateLongIDs{sIdx};
                id = substrateIDs{sIdx};
                name = substrateNames{sIdx};
                internalExchangeReactionLongIDs = [internalExchangeReactionLongIDs; {[longid 'ix']}];
                internalExchangeReactionIDs = [internalExchangeReactionIDs; {[id 'ix']}];
                internalExchangeReactionNames = [internalExchangeReactionNames; {[name ' internal exchange']}];
            end
            
            reactionLongIDs = cell(size(p.fbaReactionStoichiometryMatrix, 2), 1);
            reactionIDs = cell(size(p.fbaReactionStoichiometryMatrix, 2), 1);
            reactionNames = cell(size(p.fbaReactionStoichiometryMatrix, 2), 1);
            reversible = false(size(p.fbaReactionStoichiometryMatrix, 2), 1);
            reactionLongIDs(p.fbaReactionIndexs_metabolicConversion) = tmp_reactionLongIDs(p.reactionIndexs_fba);
            reactionIDs(p.fbaReactionIndexs_metabolicConversion) = tmp_reactionIDs(p.reactionIndexs_fba);
            reactionNames(p.fbaReactionIndexs_metabolicConversion) = tmp_reactionNames(p.reactionIndexs_fba);
            reversible(p.fbaReactionIndexs_metabolicConversion) = tmp_reversible(p.reactionIndexs_fba);
            reactionLongIDs(p.fbaReactionIndexs_metaboliteExternalExchange) = externalExchangeReactionLongIDs;
            reactionIDs(p.fbaReactionIndexs_metaboliteExternalExchange) = externalExchangeReactionIDs;
            reactionNames(p.fbaReactionIndexs_metaboliteExternalExchange) = externalExchangeReactionNames;
            reactionLongIDs(p.fbaReactionIndexs_metaboliteInternalExchange) = internalExchangeReactionLongIDs;
            reactionIDs(p.fbaReactionIndexs_metaboliteInternalExchange) = internalExchangeReactionIDs;
            reactionNames(p.fbaReactionIndexs_metaboliteInternalExchange) = internalExchangeReactionNames;
            reactionLongIDs{p.fbaReactionIndexs_biomassProduction} = 'Growth';
            reactionIDs{p.fbaReactionIndexs_biomassProduction} = 'Growth';
            reactionNames{p.fbaReactionIndexs_biomassProduction} = 'Growth';
            reactionLongIDs{p.fbaReactionIndexs_biomassExchange} = 'Bmssex';
            reactionIDs{p.fbaReactionIndexs_biomassExchange} = 'Bmssex';
            reactionNames{p.fbaReactionIndexs_biomassExchange} = 'Biomass exchange';
            
            assertEqual(reversible, upper(cellfun(@(id) id(end), reactionIDs)) == 'R');
            assertEqual(numel(reactionIDs), numel(unique(reactionIDs)));
            assertIn(max(cellfun(@length, reactionIDs)), [0 8]);
            
            %% substrates
            if formatSubstrateIdsForSbml
                substrateLongIDs(~cellfun(@isempty, substrateCompartmentIDs)) = ...
                    cellfun(@(x,y) [x '[' y ']'], ...
                    substrateLongIDs(~cellfun(@isempty, substrateCompartmentIDs)), ...
                    substrateCompartmentIDs(~cellfun(@isempty, substrateCompartmentIDs)), ...
                    'UniformOutput', false);
                substrateLongIDs(cellfun(@isempty, substrateCompartmentIDs)) = ...
                    cellfun(@(x) [x '[none]'], ...
                    substrateLongIDs(cellfun(@isempty, substrateCompartmentIDs)), ...
                    'UniformOutput', false);
            else
                substrateLongIDs(~cellfun(@isempty, substrateCompartmentIDs)) = ...
                    cellfun(@(x,y) [x '_' y], ...
                    substrateLongIDs(~cellfun(@isempty, substrateCompartmentIDs)), ...
                    substrateCompartmentIDs(~cellfun(@isempty, substrateCompartmentIDs)), ...
                    'UniformOutput', false);
            end
            substrateIDs = cellfun(@(x,y) [x y], substrateIDs, substrateCompartmentIDs, 'UniformOutput', false);
        end
    end
end
