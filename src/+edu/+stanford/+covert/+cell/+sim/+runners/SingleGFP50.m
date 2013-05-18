classdef SingleGFP50 < edu.stanford.covert.cell.sim.runners.SimulationRunner
	methods
		function this = SingleGFP50(varargin)
			this = this@edu.stanford.covert.cell.sim.runners.SimulationRunner(varargin{:});
		end
	end

	methods (Access = protected)
		function modifyNetworkStructure(this, kb)
			import edu.stanford.covert.cell.kb.Gene;
			import edu.stanford.covert.cell.kb.ProteinMonomer;
			import edu.stanford.covert.cell.kb.ProteinComplex;
			import edu.stanford.covert.cell.kb.Stimuli;
			import edu.stanford.covert.cell.kb.TranscriptionUnit;
			meanHL = mean([kb.mRNAGenes.halfLife]);
			expression = zeros(1, 3);
			geneA1 = Gene(kb, NaN, {'GeneA1'}, {'Gene A1'}, {'genA1'}, {''}, {'mRNA'}, 0, {''}, 580177, 717, true, {''}, meanHL, expression);
			geneA2 = Gene(kb, NaN, {'GeneA2'}, {'Gene A2'}, {'genA2'}, {''}, {'mRNA'}, 0, {''}, 580994, 717, true, {''}, meanHL, expression);
			geneA3 = Gene(kb, NaN, {'GeneA3'}, {'Gene A3'}, {'genA3'}, {''}, {'mRNA'}, 0, {''}, 581811, 717, true, {''}, meanHL, expression);
			geneA4 = Gene(kb, NaN, {'GeneA4'}, {'Gene A4'}, {'genA4'}, {''}, {'mRNA'}, 0, {''}, 582628, 717, true, {''}, meanHL, expression);
			geneA5 = Gene(kb, NaN, {'GeneA5'}, {'Gene A5'}, {'genA5'}, {''}, {'mRNA'}, 0, {''}, 583445, 717, true, {''}, meanHL, expression);
			geneA6 = Gene(kb, NaN, {'GeneA6'}, {'Gene A6'}, {'genA6'}, {''}, {'mRNA'}, 0, {''}, 584262, 717, true, {''}, meanHL, expression);
			geneA7 = Gene(kb, NaN, {'GeneA7'}, {'Gene A7'}, {'genA7'}, {''}, {'mRNA'}, 0, {''}, 585079, 717, true, {''}, meanHL, expression);
			geneA8 = Gene(kb, NaN, {'GeneA8'}, {'Gene A8'}, {'genA8'}, {''}, {'mRNA'}, 0, {''}, 585896, 717, true, {''}, meanHL, expression);
			geneA9 = Gene(kb, NaN, {'GeneA9'}, {'Gene A9'}, {'genA9'}, {''}, {'mRNA'}, 0, {''}, 586713, 717, true, {''}, meanHL, expression);
			geneA10 = Gene(kb, NaN, {'GeneA10'}, {'Gene A10'}, {'genA10'}, {''}, {'mRNA'}, 0, {''}, 587530, 717, true, {''}, meanHL, expression);
			geneA11 = Gene(kb, NaN, {'GeneA11'}, {'Gene A11'}, {'genA11'}, {''}, {'mRNA'}, 0, {''}, 588347, 717, true, {''}, meanHL, expression);
			geneA12 = Gene(kb, NaN, {'GeneA12'}, {'Gene A12'}, {'genA12'}, {''}, {'mRNA'}, 0, {''}, 589164, 717, true, {''}, meanHL, expression);
			geneA13 = Gene(kb, NaN, {'GeneA13'}, {'Gene A13'}, {'genA13'}, {''}, {'mRNA'}, 0, {''}, 589981, 717, true, {''}, meanHL, expression);
			geneA14 = Gene(kb, NaN, {'GeneA14'}, {'Gene A14'}, {'genA14'}, {''}, {'mRNA'}, 0, {''}, 590798, 717, true, {''}, meanHL, expression);
			geneA15 = Gene(kb, NaN, {'GeneA15'}, {'Gene A15'}, {'genA15'}, {''}, {'mRNA'}, 0, {''}, 591615, 717, true, {''}, meanHL, expression);
			geneA16 = Gene(kb, NaN, {'GeneA16'}, {'Gene A16'}, {'genA16'}, {''}, {'mRNA'}, 0, {''}, 592432, 717, true, {''}, meanHL, expression);
			geneA17 = Gene(kb, NaN, {'GeneA17'}, {'Gene A17'}, {'genA17'}, {''}, {'mRNA'}, 0, {''}, 593249, 717, true, {''}, meanHL, expression);
			geneA18 = Gene(kb, NaN, {'GeneA18'}, {'Gene A18'}, {'genA18'}, {''}, {'mRNA'}, 0, {''}, 594066, 717, true, {''}, meanHL, expression);
			geneA19 = Gene(kb, NaN, {'GeneA19'}, {'Gene A19'}, {'genA19'}, {''}, {'mRNA'}, 0, {''}, 594883, 717, true, {''}, meanHL, expression);
			geneA20 = Gene(kb, NaN, {'GeneA20'}, {'Gene A20'}, {'genA20'}, {''}, {'mRNA'}, 0, {''}, 595700, 717, true, {''}, meanHL, expression);
			geneA21 = Gene(kb, NaN, {'GeneA21'}, {'Gene A21'}, {'genA21'}, {''}, {'mRNA'}, 0, {''}, 596517, 717, true, {''}, meanHL, expression);
			geneA22 = Gene(kb, NaN, {'GeneA22'}, {'Gene A22'}, {'genA22'}, {''}, {'mRNA'}, 0, {''}, 597334, 717, true, {''}, meanHL, expression);
			geneA23 = Gene(kb, NaN, {'GeneA23'}, {'Gene A23'}, {'genA23'}, {''}, {'mRNA'}, 0, {''}, 598151, 717, true, {''}, meanHL, expression);
			geneA24 = Gene(kb, NaN, {'GeneA24'}, {'Gene A24'}, {'genA24'}, {''}, {'mRNA'}, 0, {''}, 598968, 717, true, {''}, meanHL, expression);
			geneA25 = Gene(kb, NaN, {'GeneA25'}, {'Gene A25'}, {'genA25'}, {''}, {'mRNA'}, 0, {''}, 599785, 717, true, {''}, meanHL, expression);
			geneA26 = Gene(kb, NaN, {'GeneA26'}, {'Gene A26'}, {'genA26'}, {''}, {'mRNA'}, 0, {''}, 600602, 717, true, {''}, meanHL, expression);
			geneA27 = Gene(kb, NaN, {'GeneA27'}, {'Gene A27'}, {'genA27'}, {''}, {'mRNA'}, 0, {''}, 601419, 717, true, {''}, meanHL, expression);
			geneA28 = Gene(kb, NaN, {'GeneA28'}, {'Gene A28'}, {'genA28'}, {''}, {'mRNA'}, 0, {''}, 602236, 717, true, {''}, meanHL, expression);
			geneA29 = Gene(kb, NaN, {'GeneA29'}, {'Gene A29'}, {'genA29'}, {''}, {'mRNA'}, 0, {''}, 603053, 717, true, {''}, meanHL, expression);
			geneA30 = Gene(kb, NaN, {'GeneA30'}, {'Gene A30'}, {'genA30'}, {''}, {'mRNA'}, 0, {''}, 603870, 717, true, {''}, meanHL, expression);
			geneA31 = Gene(kb, NaN, {'GeneA31'}, {'Gene A31'}, {'genA31'}, {''}, {'mRNA'}, 0, {''}, 604687, 717, true, {''}, meanHL, expression);
			geneA32 = Gene(kb, NaN, {'GeneA32'}, {'Gene A32'}, {'genA32'}, {''}, {'mRNA'}, 0, {''}, 605504, 717, true, {''}, meanHL, expression);
			geneA33 = Gene(kb, NaN, {'GeneA33'}, {'Gene A33'}, {'genA33'}, {''}, {'mRNA'}, 0, {''}, 606321, 717, true, {''}, meanHL, expression);
			geneA34 = Gene(kb, NaN, {'GeneA34'}, {'Gene A34'}, {'genA34'}, {''}, {'mRNA'}, 0, {''}, 607138, 717, true, {''}, meanHL, expression);
			geneA35 = Gene(kb, NaN, {'GeneA35'}, {'Gene A35'}, {'genA35'}, {''}, {'mRNA'}, 0, {''}, 607955, 717, true, {''}, meanHL, expression);
			geneA36 = Gene(kb, NaN, {'GeneA36'}, {'Gene A36'}, {'genA36'}, {''}, {'mRNA'}, 0, {''}, 608772, 717, true, {''}, meanHL, expression);
			geneA37 = Gene(kb, NaN, {'GeneA37'}, {'Gene A37'}, {'genA37'}, {''}, {'mRNA'}, 0, {''}, 609589, 717, true, {''}, meanHL, expression);
			geneA38 = Gene(kb, NaN, {'GeneA38'}, {'Gene A38'}, {'genA38'}, {''}, {'mRNA'}, 0, {''}, 610406, 717, true, {''}, meanHL, expression);
			geneA39 = Gene(kb, NaN, {'GeneA39'}, {'Gene A39'}, {'genA39'}, {''}, {'mRNA'}, 0, {''}, 611223, 717, true, {''}, meanHL, expression);
			geneA40 = Gene(kb, NaN, {'GeneA40'}, {'Gene A40'}, {'genA40'}, {''}, {'mRNA'}, 0, {''}, 612040, 717, true, {''}, meanHL, expression);
			geneA41 = Gene(kb, NaN, {'GeneA41'}, {'Gene A41'}, {'genA41'}, {''}, {'mRNA'}, 0, {''}, 612857, 717, true, {''}, meanHL, expression);
			geneA42 = Gene(kb, NaN, {'GeneA42'}, {'Gene A42'}, {'genA42'}, {''}, {'mRNA'}, 0, {''}, 613674, 717, true, {''}, meanHL, expression);
			geneA43 = Gene(kb, NaN, {'GeneA43'}, {'Gene A43'}, {'genA43'}, {''}, {'mRNA'}, 0, {''}, 614491, 717, true, {''}, meanHL, expression);
			geneA44 = Gene(kb, NaN, {'GeneA44'}, {'Gene A44'}, {'genA44'}, {''}, {'mRNA'}, 0, {''}, 615308, 717, true, {''}, meanHL, expression);
			geneA45 = Gene(kb, NaN, {'GeneA45'}, {'Gene A45'}, {'genA45'}, {''}, {'mRNA'}, 0, {''}, 616125, 717, true, {''}, meanHL, expression);
			geneA46 = Gene(kb, NaN, {'GeneA46'}, {'Gene A46'}, {'genA46'}, {''}, {'mRNA'}, 0, {''}, 616942, 717, true, {''}, meanHL, expression);
			geneA47 = Gene(kb, NaN, {'GeneA47'}, {'Gene A47'}, {'genA47'}, {''}, {'mRNA'}, 0, {''}, 617759, 717, true, {''}, meanHL, expression);
			geneA48 = Gene(kb, NaN, {'GeneA48'}, {'Gene A48'}, {'genA48'}, {''}, {'mRNA'}, 0, {''}, 618576, 717, true, {''}, meanHL, expression);
			geneA49 = Gene(kb, NaN, {'GeneA49'}, {'Gene A49'}, {'genA49'}, {''}, {'mRNA'}, 0, {''}, 619393, 717, true, {''}, meanHL, expression);
			geneA50 = Gene(kb, NaN, {'GeneA50'}, {'Gene A50'}, {'genA50'}, {''}, {'mRNA'}, 0, {''}, 620210, 717, true, {''}, meanHL, expression);
			tuA1 = TranscriptionUnit(kb, NaN, {'TuA1'}, {'Transcription unit A1'}, -35, 6, -35, 6, -1);
			tuA2 = TranscriptionUnit(kb, NaN, {'TuA2'}, {'Transcription unit A2'}, -35, 6, -35, 6, -1);
			tuA3 = TranscriptionUnit(kb, NaN, {'TuA3'}, {'Transcription unit A3'}, -35, 6, -35, 6, -1);
			tuA4 = TranscriptionUnit(kb, NaN, {'TuA4'}, {'Transcription unit A4'}, -35, 6, -35, 6, -1);
			tuA5 = TranscriptionUnit(kb, NaN, {'TuA5'}, {'Transcription unit A5'}, -35, 6, -35, 6, -1);
			tuA6 = TranscriptionUnit(kb, NaN, {'TuA6'}, {'Transcription unit A6'}, -35, 6, -35, 6, -1);
			tuA7 = TranscriptionUnit(kb, NaN, {'TuA7'}, {'Transcription unit A7'}, -35, 6, -35, 6, -1);
			tuA8 = TranscriptionUnit(kb, NaN, {'TuA8'}, {'Transcription unit A8'}, -35, 6, -35, 6, -1);
			tuA9 = TranscriptionUnit(kb, NaN, {'TuA9'}, {'Transcription unit A9'}, -35, 6, -35, 6, -1);
			tuA10 = TranscriptionUnit(kb, NaN, {'TuA10'}, {'Transcription unit A10'}, -35, 6, -35, 6, -1);
			tuA11 = TranscriptionUnit(kb, NaN, {'TuA11'}, {'Transcription unit A11'}, -35, 6, -35, 6, -1);
			tuA12 = TranscriptionUnit(kb, NaN, {'TuA12'}, {'Transcription unit A12'}, -35, 6, -35, 6, -1);
			tuA13 = TranscriptionUnit(kb, NaN, {'TuA13'}, {'Transcription unit A13'}, -35, 6, -35, 6, -1);
			tuA14 = TranscriptionUnit(kb, NaN, {'TuA14'}, {'Transcription unit A14'}, -35, 6, -35, 6, -1);
			tuA15 = TranscriptionUnit(kb, NaN, {'TuA15'}, {'Transcription unit A15'}, -35, 6, -35, 6, -1);
			tuA16 = TranscriptionUnit(kb, NaN, {'TuA16'}, {'Transcription unit A16'}, -35, 6, -35, 6, -1);
			tuA17 = TranscriptionUnit(kb, NaN, {'TuA17'}, {'Transcription unit A17'}, -35, 6, -35, 6, -1);
			tuA18 = TranscriptionUnit(kb, NaN, {'TuA18'}, {'Transcription unit A18'}, -35, 6, -35, 6, -1);
			tuA19 = TranscriptionUnit(kb, NaN, {'TuA19'}, {'Transcription unit A19'}, -35, 6, -35, 6, -1);
			tuA20 = TranscriptionUnit(kb, NaN, {'TuA20'}, {'Transcription unit A20'}, -35, 6, -35, 6, -1);
			tuA21 = TranscriptionUnit(kb, NaN, {'TuA21'}, {'Transcription unit A21'}, -35, 6, -35, 6, -1);
			tuA22 = TranscriptionUnit(kb, NaN, {'TuA22'}, {'Transcription unit A22'}, -35, 6, -35, 6, -1);
			tuA23 = TranscriptionUnit(kb, NaN, {'TuA23'}, {'Transcription unit A23'}, -35, 6, -35, 6, -1);
			tuA24 = TranscriptionUnit(kb, NaN, {'TuA24'}, {'Transcription unit A24'}, -35, 6, -35, 6, -1);
			tuA25 = TranscriptionUnit(kb, NaN, {'TuA25'}, {'Transcription unit A25'}, -35, 6, -35, 6, -1);
			tuA26 = TranscriptionUnit(kb, NaN, {'TuA26'}, {'Transcription unit A26'}, -35, 6, -35, 6, -1);
			tuA27 = TranscriptionUnit(kb, NaN, {'TuA27'}, {'Transcription unit A27'}, -35, 6, -35, 6, -1);
			tuA28 = TranscriptionUnit(kb, NaN, {'TuA28'}, {'Transcription unit A28'}, -35, 6, -35, 6, -1);
			tuA29 = TranscriptionUnit(kb, NaN, {'TuA29'}, {'Transcription unit A29'}, -35, 6, -35, 6, -1);
			tuA30 = TranscriptionUnit(kb, NaN, {'TuA30'}, {'Transcription unit A30'}, -35, 6, -35, 6, -1);
			tuA31 = TranscriptionUnit(kb, NaN, {'TuA31'}, {'Transcription unit A31'}, -35, 6, -35, 6, -1);
			tuA32 = TranscriptionUnit(kb, NaN, {'TuA32'}, {'Transcription unit A32'}, -35, 6, -35, 6, -1);
			tuA33 = TranscriptionUnit(kb, NaN, {'TuA33'}, {'Transcription unit A33'}, -35, 6, -35, 6, -1);
			tuA34 = TranscriptionUnit(kb, NaN, {'TuA34'}, {'Transcription unit A34'}, -35, 6, -35, 6, -1);
			tuA35 = TranscriptionUnit(kb, NaN, {'TuA35'}, {'Transcription unit A35'}, -35, 6, -35, 6, -1);
			tuA36 = TranscriptionUnit(kb, NaN, {'TuA36'}, {'Transcription unit A36'}, -35, 6, -35, 6, -1);
			tuA37 = TranscriptionUnit(kb, NaN, {'TuA37'}, {'Transcription unit A37'}, -35, 6, -35, 6, -1);
			tuA38 = TranscriptionUnit(kb, NaN, {'TuA38'}, {'Transcription unit A38'}, -35, 6, -35, 6, -1);
			tuA39 = TranscriptionUnit(kb, NaN, {'TuA39'}, {'Transcription unit A39'}, -35, 6, -35, 6, -1);
			tuA40 = TranscriptionUnit(kb, NaN, {'TuA40'}, {'Transcription unit A40'}, -35, 6, -35, 6, -1);
			tuA41 = TranscriptionUnit(kb, NaN, {'TuA41'}, {'Transcription unit A41'}, -35, 6, -35, 6, -1);
			tuA42 = TranscriptionUnit(kb, NaN, {'TuA42'}, {'Transcription unit A42'}, -35, 6, -35, 6, -1);
			tuA43 = TranscriptionUnit(kb, NaN, {'TuA43'}, {'Transcription unit A43'}, -35, 6, -35, 6, -1);
			tuA44 = TranscriptionUnit(kb, NaN, {'TuA44'}, {'Transcription unit A44'}, -35, 6, -35, 6, -1);
			tuA45 = TranscriptionUnit(kb, NaN, {'TuA45'}, {'Transcription unit A45'}, -35, 6, -35, 6, -1);
			tuA46 = TranscriptionUnit(kb, NaN, {'TuA46'}, {'Transcription unit A46'}, -35, 6, -35, 6, -1);
			tuA47 = TranscriptionUnit(kb, NaN, {'TuA47'}, {'Transcription unit A47'}, -35, 6, -35, 6, -1);
			tuA48 = TranscriptionUnit(kb, NaN, {'TuA48'}, {'Transcription unit A48'}, -35, 6, -35, 6, -1);
			tuA49 = TranscriptionUnit(kb, NaN, {'TuA49'}, {'Transcription unit A49'}, -35, 6, -35, 6, -1);
			tuA50 = TranscriptionUnit(kb, NaN, {'TuA50'}, {'Transcription unit A50'}, -35, 6, -35, 6, -1);
			monA1 = ProteinMonomer(kb, NaN, {'MonomerA1'}, {'Protein monomer A1'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA2 = ProteinMonomer(kb, NaN, {'MonomerA2'}, {'Protein monomer A2'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA3 = ProteinMonomer(kb, NaN, {'MonomerA3'}, {'Protein monomer A3'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA4 = ProteinMonomer(kb, NaN, {'MonomerA4'}, {'Protein monomer A4'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA5 = ProteinMonomer(kb, NaN, {'MonomerA5'}, {'Protein monomer A5'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA6 = ProteinMonomer(kb, NaN, {'MonomerA6'}, {'Protein monomer A6'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA7 = ProteinMonomer(kb, NaN, {'MonomerA7'}, {'Protein monomer A7'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA8 = ProteinMonomer(kb, NaN, {'MonomerA8'}, {'Protein monomer A8'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA9 = ProteinMonomer(kb, NaN, {'MonomerA9'}, {'Protein monomer A9'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA10 = ProteinMonomer(kb, NaN, {'MonomerA10'}, {'Protein monomer A10'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA11 = ProteinMonomer(kb, NaN, {'MonomerA11'}, {'Protein monomer A11'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA12 = ProteinMonomer(kb, NaN, {'MonomerA12'}, {'Protein monomer A12'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA13 = ProteinMonomer(kb, NaN, {'MonomerA13'}, {'Protein monomer A13'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA14 = ProteinMonomer(kb, NaN, {'MonomerA14'}, {'Protein monomer A14'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA15 = ProteinMonomer(kb, NaN, {'MonomerA15'}, {'Protein monomer A15'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA16 = ProteinMonomer(kb, NaN, {'MonomerA16'}, {'Protein monomer A16'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA17 = ProteinMonomer(kb, NaN, {'MonomerA17'}, {'Protein monomer A17'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA18 = ProteinMonomer(kb, NaN, {'MonomerA18'}, {'Protein monomer A18'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA19 = ProteinMonomer(kb, NaN, {'MonomerA19'}, {'Protein monomer A19'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA20 = ProteinMonomer(kb, NaN, {'MonomerA20'}, {'Protein monomer A20'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA21 = ProteinMonomer(kb, NaN, {'MonomerA21'}, {'Protein monomer A21'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA22 = ProteinMonomer(kb, NaN, {'MonomerA22'}, {'Protein monomer A22'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA23 = ProteinMonomer(kb, NaN, {'MonomerA23'}, {'Protein monomer A23'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA24 = ProteinMonomer(kb, NaN, {'MonomerA24'}, {'Protein monomer A24'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA25 = ProteinMonomer(kb, NaN, {'MonomerA25'}, {'Protein monomer A25'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA26 = ProteinMonomer(kb, NaN, {'MonomerA26'}, {'Protein monomer A26'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA27 = ProteinMonomer(kb, NaN, {'MonomerA27'}, {'Protein monomer A27'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA28 = ProteinMonomer(kb, NaN, {'MonomerA28'}, {'Protein monomer A28'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA29 = ProteinMonomer(kb, NaN, {'MonomerA29'}, {'Protein monomer A29'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA30 = ProteinMonomer(kb, NaN, {'MonomerA30'}, {'Protein monomer A30'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA31 = ProteinMonomer(kb, NaN, {'MonomerA31'}, {'Protein monomer A31'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA32 = ProteinMonomer(kb, NaN, {'MonomerA32'}, {'Protein monomer A32'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA33 = ProteinMonomer(kb, NaN, {'MonomerA33'}, {'Protein monomer A33'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA34 = ProteinMonomer(kb, NaN, {'MonomerA34'}, {'Protein monomer A34'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA35 = ProteinMonomer(kb, NaN, {'MonomerA35'}, {'Protein monomer A35'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA36 = ProteinMonomer(kb, NaN, {'MonomerA36'}, {'Protein monomer A36'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA37 = ProteinMonomer(kb, NaN, {'MonomerA37'}, {'Protein monomer A37'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA38 = ProteinMonomer(kb, NaN, {'MonomerA38'}, {'Protein monomer A38'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA39 = ProteinMonomer(kb, NaN, {'MonomerA39'}, {'Protein monomer A39'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA40 = ProteinMonomer(kb, NaN, {'MonomerA40'}, {'Protein monomer A40'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA41 = ProteinMonomer(kb, NaN, {'MonomerA41'}, {'Protein monomer A41'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA42 = ProteinMonomer(kb, NaN, {'MonomerA42'}, {'Protein monomer A42'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA43 = ProteinMonomer(kb, NaN, {'MonomerA43'}, {'Protein monomer A43'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA44 = ProteinMonomer(kb, NaN, {'MonomerA44'}, {'Protein monomer A44'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA45 = ProteinMonomer(kb, NaN, {'MonomerA45'}, {'Protein monomer A45'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA46 = ProteinMonomer(kb, NaN, {'MonomerA46'}, {'Protein monomer A46'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA47 = ProteinMonomer(kb, NaN, {'MonomerA47'}, {'Protein monomer A47'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA48 = ProteinMonomer(kb, NaN, {'MonomerA48'}, {'Protein monomer A48'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA49 = ProteinMonomer(kb, NaN, {'MonomerA49'}, {'Protein monomer A49'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA50 = ProteinMonomer(kb, NaN, {'MonomerA50'}, {'Protein monomer A50'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			cpxAA1 = ProteinComplex(kb, NaN, {'ComplexAA1'}, {'Macromolecular complex AA1'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA2 = ProteinComplex(kb, NaN, {'ComplexAA2'}, {'Macromolecular complex AA2'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA3 = ProteinComplex(kb, NaN, {'ComplexAA3'}, {'Macromolecular complex AA3'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA4 = ProteinComplex(kb, NaN, {'ComplexAA4'}, {'Macromolecular complex AA4'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA5 = ProteinComplex(kb, NaN, {'ComplexAA5'}, {'Macromolecular complex AA5'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA6 = ProteinComplex(kb, NaN, {'ComplexAA6'}, {'Macromolecular complex AA6'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA7 = ProteinComplex(kb, NaN, {'ComplexAA7'}, {'Macromolecular complex AA7'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA8 = ProteinComplex(kb, NaN, {'ComplexAA8'}, {'Macromolecular complex AA8'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA9 = ProteinComplex(kb, NaN, {'ComplexAA9'}, {'Macromolecular complex AA9'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA10 = ProteinComplex(kb, NaN, {'ComplexAA10'}, {'Macromolecular complex AA10'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA11 = ProteinComplex(kb, NaN, {'ComplexAA11'}, {'Macromolecular complex AA11'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA12 = ProteinComplex(kb, NaN, {'ComplexAA12'}, {'Macromolecular complex AA12'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA13 = ProteinComplex(kb, NaN, {'ComplexAA13'}, {'Macromolecular complex AA13'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA14 = ProteinComplex(kb, NaN, {'ComplexAA14'}, {'Macromolecular complex AA14'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA15 = ProteinComplex(kb, NaN, {'ComplexAA15'}, {'Macromolecular complex AA15'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA16 = ProteinComplex(kb, NaN, {'ComplexAA16'}, {'Macromolecular complex AA16'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA17 = ProteinComplex(kb, NaN, {'ComplexAA17'}, {'Macromolecular complex AA17'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA18 = ProteinComplex(kb, NaN, {'ComplexAA18'}, {'Macromolecular complex AA18'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA19 = ProteinComplex(kb, NaN, {'ComplexAA19'}, {'Macromolecular complex AA19'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA20 = ProteinComplex(kb, NaN, {'ComplexAA20'}, {'Macromolecular complex AA20'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA21 = ProteinComplex(kb, NaN, {'ComplexAA21'}, {'Macromolecular complex AA21'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA22 = ProteinComplex(kb, NaN, {'ComplexAA22'}, {'Macromolecular complex AA22'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA23 = ProteinComplex(kb, NaN, {'ComplexAA23'}, {'Macromolecular complex AA23'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA24 = ProteinComplex(kb, NaN, {'ComplexAA24'}, {'Macromolecular complex AA24'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA25 = ProteinComplex(kb, NaN, {'ComplexAA25'}, {'Macromolecular complex AA25'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA26 = ProteinComplex(kb, NaN, {'ComplexAA26'}, {'Macromolecular complex AA26'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA27 = ProteinComplex(kb, NaN, {'ComplexAA27'}, {'Macromolecular complex AA27'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA28 = ProteinComplex(kb, NaN, {'ComplexAA28'}, {'Macromolecular complex AA28'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA29 = ProteinComplex(kb, NaN, {'ComplexAA29'}, {'Macromolecular complex AA29'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA30 = ProteinComplex(kb, NaN, {'ComplexAA30'}, {'Macromolecular complex AA30'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA31 = ProteinComplex(kb, NaN, {'ComplexAA31'}, {'Macromolecular complex AA31'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA32 = ProteinComplex(kb, NaN, {'ComplexAA32'}, {'Macromolecular complex AA32'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA33 = ProteinComplex(kb, NaN, {'ComplexAA33'}, {'Macromolecular complex AA33'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA34 = ProteinComplex(kb, NaN, {'ComplexAA34'}, {'Macromolecular complex AA34'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA35 = ProteinComplex(kb, NaN, {'ComplexAA35'}, {'Macromolecular complex AA35'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA36 = ProteinComplex(kb, NaN, {'ComplexAA36'}, {'Macromolecular complex AA36'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA37 = ProteinComplex(kb, NaN, {'ComplexAA37'}, {'Macromolecular complex AA37'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA38 = ProteinComplex(kb, NaN, {'ComplexAA38'}, {'Macromolecular complex AA38'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA39 = ProteinComplex(kb, NaN, {'ComplexAA39'}, {'Macromolecular complex AA39'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA40 = ProteinComplex(kb, NaN, {'ComplexAA40'}, {'Macromolecular complex AA40'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA41 = ProteinComplex(kb, NaN, {'ComplexAA41'}, {'Macromolecular complex AA41'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA42 = ProteinComplex(kb, NaN, {'ComplexAA42'}, {'Macromolecular complex AA42'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA43 = ProteinComplex(kb, NaN, {'ComplexAA43'}, {'Macromolecular complex AA43'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA44 = ProteinComplex(kb, NaN, {'ComplexAA44'}, {'Macromolecular complex AA44'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA45 = ProteinComplex(kb, NaN, {'ComplexAA45'}, {'Macromolecular complex AA45'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA46 = ProteinComplex(kb, NaN, {'ComplexAA46'}, {'Macromolecular complex AA46'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA47 = ProteinComplex(kb, NaN, {'ComplexAA47'}, {'Macromolecular complex AA47'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA48 = ProteinComplex(kb, NaN, {'ComplexAA48'}, {'Macromolecular complex AA48'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA49 = ProteinComplex(kb, NaN, {'ComplexAA49'}, {'Macromolecular complex AA49'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA50 = ProteinComplex(kb, NaN, {'ComplexAA50'}, {'Macromolecular complex AA50'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			genome = kb.genome;
			promoter = 'TATCTACAGAGGCCCGCTTTAAGCATCCACGACCCTACTCACTTCAAAGTGGAACCACCCGTCGACGTGTGTCTTAACCCCTGCCGCTGCAAGGTGTGAG';
			gfp = 'ATGTCTAAAGGTGAAGAATTATTCACTGGTGTTGTCCCAATTTTGGTTGAATTAGATGGTGATGTTAATGGTCACAAATTTTCTGTCTCCGGTGAAGGTGAAGGTGATGCTACTTACGGTAAATTGACCTTAAAATTTATTTGTACTACTGGTAAATTGCCAGTTCCATGGCCAACCTTAGTCACTACTTTCGGTTATGGTGTTCAATGTTTTGCTAGATACCCAGATCATATGAAACAACATGACTTTTTCAAGTCTGCCATGCCAGAAGGTTATGTTCAAGAAAGAACTATTTTTTTCAAAGATGACGGTAACTACAAGACCAGAGCTGAAGTCAAGTTTGAAGGTGATACCTTAGTTAATAGAATCGAATTAAAAGGTATTGATTTTAAAGAAGATGGTAACATTTTAGGTCACAAATTGGAATACAACTATAACTCTCACAATGTTTACATCATGGCTGACAAACAAAAGAATGGTATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGTTCTGTTCAATTAGCTGACCATTATCAACAAAATACTCCAATTGGTGATGGTCCAGTCTTGTTACCAGACAACCATTACTTATCCACTCAATCTGCCTTATCCAAAGATCCAAACGAAAAGAGAGACCACATGGTCTTGTTAGAATTTGTTACTGCTGCTGGTATTACCCATGGTATGGATGAATTGTACAAATAA';
			genome.sequence = [genome.sequence repmat([promoter gfp], [1, 50])];
			genome.sequenceLength = length(genome.sequence);
			cytosol = kb.compartments(kb.cytosolCompartmentIndexs);
			extracellularSpace = kb.compartments(kb.extracellularCompartmentIndexs);
			geneA1.genome = genome;
			geneA1.transcriptionUnits = tuA1;
			geneA1.proteinMonomers = monA1;
			geneA1.compartment = cytosol;

			geneA2.genome = genome;
			geneA2.transcriptionUnits = tuA2;
			geneA2.proteinMonomers = monA2;
			geneA2.compartment = cytosol;

			geneA3.genome = genome;
			geneA3.transcriptionUnits = tuA3;
			geneA3.proteinMonomers = monA3;
			geneA3.compartment = cytosol;

			geneA4.genome = genome;
			geneA4.transcriptionUnits = tuA4;
			geneA4.proteinMonomers = monA4;
			geneA4.compartment = cytosol;

			geneA5.genome = genome;
			geneA5.transcriptionUnits = tuA5;
			geneA5.proteinMonomers = monA5;
			geneA5.compartment = cytosol;

			geneA6.genome = genome;
			geneA6.transcriptionUnits = tuA6;
			geneA6.proteinMonomers = monA6;
			geneA6.compartment = cytosol;

			geneA7.genome = genome;
			geneA7.transcriptionUnits = tuA7;
			geneA7.proteinMonomers = monA7;
			geneA7.compartment = cytosol;

			geneA8.genome = genome;
			geneA8.transcriptionUnits = tuA8;
			geneA8.proteinMonomers = monA8;
			geneA8.compartment = cytosol;

			geneA9.genome = genome;
			geneA9.transcriptionUnits = tuA9;
			geneA9.proteinMonomers = monA9;
			geneA9.compartment = cytosol;

			geneA10.genome = genome;
			geneA10.transcriptionUnits = tuA10;
			geneA10.proteinMonomers = monA10;
			geneA10.compartment = cytosol;

			geneA11.genome = genome;
			geneA11.transcriptionUnits = tuA11;
			geneA11.proteinMonomers = monA11;
			geneA11.compartment = cytosol;

			geneA12.genome = genome;
			geneA12.transcriptionUnits = tuA12;
			geneA12.proteinMonomers = monA12;
			geneA12.compartment = cytosol;

			geneA13.genome = genome;
			geneA13.transcriptionUnits = tuA13;
			geneA13.proteinMonomers = monA13;
			geneA13.compartment = cytosol;

			geneA14.genome = genome;
			geneA14.transcriptionUnits = tuA14;
			geneA14.proteinMonomers = monA14;
			geneA14.compartment = cytosol;

			geneA15.genome = genome;
			geneA15.transcriptionUnits = tuA15;
			geneA15.proteinMonomers = monA15;
			geneA15.compartment = cytosol;

			geneA16.genome = genome;
			geneA16.transcriptionUnits = tuA16;
			geneA16.proteinMonomers = monA16;
			geneA16.compartment = cytosol;

			geneA17.genome = genome;
			geneA17.transcriptionUnits = tuA17;
			geneA17.proteinMonomers = monA17;
			geneA17.compartment = cytosol;

			geneA18.genome = genome;
			geneA18.transcriptionUnits = tuA18;
			geneA18.proteinMonomers = monA18;
			geneA18.compartment = cytosol;

			geneA19.genome = genome;
			geneA19.transcriptionUnits = tuA19;
			geneA19.proteinMonomers = monA19;
			geneA19.compartment = cytosol;

			geneA20.genome = genome;
			geneA20.transcriptionUnits = tuA20;
			geneA20.proteinMonomers = monA20;
			geneA20.compartment = cytosol;

			geneA21.genome = genome;
			geneA21.transcriptionUnits = tuA21;
			geneA21.proteinMonomers = monA21;
			geneA21.compartment = cytosol;

			geneA22.genome = genome;
			geneA22.transcriptionUnits = tuA22;
			geneA22.proteinMonomers = monA22;
			geneA22.compartment = cytosol;

			geneA23.genome = genome;
			geneA23.transcriptionUnits = tuA23;
			geneA23.proteinMonomers = monA23;
			geneA23.compartment = cytosol;

			geneA24.genome = genome;
			geneA24.transcriptionUnits = tuA24;
			geneA24.proteinMonomers = monA24;
			geneA24.compartment = cytosol;

			geneA25.genome = genome;
			geneA25.transcriptionUnits = tuA25;
			geneA25.proteinMonomers = monA25;
			geneA25.compartment = cytosol;

			geneA26.genome = genome;
			geneA26.transcriptionUnits = tuA26;
			geneA26.proteinMonomers = monA26;
			geneA26.compartment = cytosol;

			geneA27.genome = genome;
			geneA27.transcriptionUnits = tuA27;
			geneA27.proteinMonomers = monA27;
			geneA27.compartment = cytosol;

			geneA28.genome = genome;
			geneA28.transcriptionUnits = tuA28;
			geneA28.proteinMonomers = monA28;
			geneA28.compartment = cytosol;

			geneA29.genome = genome;
			geneA29.transcriptionUnits = tuA29;
			geneA29.proteinMonomers = monA29;
			geneA29.compartment = cytosol;

			geneA30.genome = genome;
			geneA30.transcriptionUnits = tuA30;
			geneA30.proteinMonomers = monA30;
			geneA30.compartment = cytosol;

			geneA31.genome = genome;
			geneA31.transcriptionUnits = tuA31;
			geneA31.proteinMonomers = monA31;
			geneA31.compartment = cytosol;

			geneA32.genome = genome;
			geneA32.transcriptionUnits = tuA32;
			geneA32.proteinMonomers = monA32;
			geneA32.compartment = cytosol;

			geneA33.genome = genome;
			geneA33.transcriptionUnits = tuA33;
			geneA33.proteinMonomers = monA33;
			geneA33.compartment = cytosol;

			geneA34.genome = genome;
			geneA34.transcriptionUnits = tuA34;
			geneA34.proteinMonomers = monA34;
			geneA34.compartment = cytosol;

			geneA35.genome = genome;
			geneA35.transcriptionUnits = tuA35;
			geneA35.proteinMonomers = monA35;
			geneA35.compartment = cytosol;

			geneA36.genome = genome;
			geneA36.transcriptionUnits = tuA36;
			geneA36.proteinMonomers = monA36;
			geneA36.compartment = cytosol;

			geneA37.genome = genome;
			geneA37.transcriptionUnits = tuA37;
			geneA37.proteinMonomers = monA37;
			geneA37.compartment = cytosol;

			geneA38.genome = genome;
			geneA38.transcriptionUnits = tuA38;
			geneA38.proteinMonomers = monA38;
			geneA38.compartment = cytosol;

			geneA39.genome = genome;
			geneA39.transcriptionUnits = tuA39;
			geneA39.proteinMonomers = monA39;
			geneA39.compartment = cytosol;

			geneA40.genome = genome;
			geneA40.transcriptionUnits = tuA40;
			geneA40.proteinMonomers = monA40;
			geneA40.compartment = cytosol;

			geneA41.genome = genome;
			geneA41.transcriptionUnits = tuA41;
			geneA41.proteinMonomers = monA41;
			geneA41.compartment = cytosol;

			geneA42.genome = genome;
			geneA42.transcriptionUnits = tuA42;
			geneA42.proteinMonomers = monA42;
			geneA42.compartment = cytosol;

			geneA43.genome = genome;
			geneA43.transcriptionUnits = tuA43;
			geneA43.proteinMonomers = monA43;
			geneA43.compartment = cytosol;

			geneA44.genome = genome;
			geneA44.transcriptionUnits = tuA44;
			geneA44.proteinMonomers = monA44;
			geneA44.compartment = cytosol;

			geneA45.genome = genome;
			geneA45.transcriptionUnits = tuA45;
			geneA45.proteinMonomers = monA45;
			geneA45.compartment = cytosol;

			geneA46.genome = genome;
			geneA46.transcriptionUnits = tuA46;
			geneA46.proteinMonomers = monA46;
			geneA46.compartment = cytosol;

			geneA47.genome = genome;
			geneA47.transcriptionUnits = tuA47;
			geneA47.proteinMonomers = monA47;
			geneA47.compartment = cytosol;

			geneA48.genome = genome;
			geneA48.transcriptionUnits = tuA48;
			geneA48.proteinMonomers = monA48;
			geneA48.compartment = cytosol;

			geneA49.genome = genome;
			geneA49.transcriptionUnits = tuA49;
			geneA49.proteinMonomers = monA49;
			geneA49.compartment = cytosol;

			geneA50.genome = genome;
			geneA50.transcriptionUnits = tuA50;
			geneA50.proteinMonomers = monA50;
			geneA50.compartment = cytosol;

			tuA1.genome = genome;
			tuA1.genes = geneA1;
			tuA1.geneCompartments = cytosol;
			tuA1.compartment = cytosol;

			tuA2.genome = genome;
			tuA2.genes = geneA2;
			tuA2.geneCompartments = cytosol;
			tuA2.compartment = cytosol;

			tuA3.genome = genome;
			tuA3.genes = geneA3;
			tuA3.geneCompartments = cytosol;
			tuA3.compartment = cytosol;

			tuA4.genome = genome;
			tuA4.genes = geneA4;
			tuA4.geneCompartments = cytosol;
			tuA4.compartment = cytosol;

			tuA5.genome = genome;
			tuA5.genes = geneA5;
			tuA5.geneCompartments = cytosol;
			tuA5.compartment = cytosol;

			tuA6.genome = genome;
			tuA6.genes = geneA6;
			tuA6.geneCompartments = cytosol;
			tuA6.compartment = cytosol;

			tuA7.genome = genome;
			tuA7.genes = geneA7;
			tuA7.geneCompartments = cytosol;
			tuA7.compartment = cytosol;

			tuA8.genome = genome;
			tuA8.genes = geneA8;
			tuA8.geneCompartments = cytosol;
			tuA8.compartment = cytosol;

			tuA9.genome = genome;
			tuA9.genes = geneA9;
			tuA9.geneCompartments = cytosol;
			tuA9.compartment = cytosol;

			tuA10.genome = genome;
			tuA10.genes = geneA10;
			tuA10.geneCompartments = cytosol;
			tuA10.compartment = cytosol;

			tuA11.genome = genome;
			tuA11.genes = geneA11;
			tuA11.geneCompartments = cytosol;
			tuA11.compartment = cytosol;

			tuA12.genome = genome;
			tuA12.genes = geneA12;
			tuA12.geneCompartments = cytosol;
			tuA12.compartment = cytosol;

			tuA13.genome = genome;
			tuA13.genes = geneA13;
			tuA13.geneCompartments = cytosol;
			tuA13.compartment = cytosol;

			tuA14.genome = genome;
			tuA14.genes = geneA14;
			tuA14.geneCompartments = cytosol;
			tuA14.compartment = cytosol;

			tuA15.genome = genome;
			tuA15.genes = geneA15;
			tuA15.geneCompartments = cytosol;
			tuA15.compartment = cytosol;

			tuA16.genome = genome;
			tuA16.genes = geneA16;
			tuA16.geneCompartments = cytosol;
			tuA16.compartment = cytosol;

			tuA17.genome = genome;
			tuA17.genes = geneA17;
			tuA17.geneCompartments = cytosol;
			tuA17.compartment = cytosol;

			tuA18.genome = genome;
			tuA18.genes = geneA18;
			tuA18.geneCompartments = cytosol;
			tuA18.compartment = cytosol;

			tuA19.genome = genome;
			tuA19.genes = geneA19;
			tuA19.geneCompartments = cytosol;
			tuA19.compartment = cytosol;

			tuA20.genome = genome;
			tuA20.genes = geneA20;
			tuA20.geneCompartments = cytosol;
			tuA20.compartment = cytosol;

			tuA21.genome = genome;
			tuA21.genes = geneA21;
			tuA21.geneCompartments = cytosol;
			tuA21.compartment = cytosol;

			tuA22.genome = genome;
			tuA22.genes = geneA22;
			tuA22.geneCompartments = cytosol;
			tuA22.compartment = cytosol;

			tuA23.genome = genome;
			tuA23.genes = geneA23;
			tuA23.geneCompartments = cytosol;
			tuA23.compartment = cytosol;

			tuA24.genome = genome;
			tuA24.genes = geneA24;
			tuA24.geneCompartments = cytosol;
			tuA24.compartment = cytosol;

			tuA25.genome = genome;
			tuA25.genes = geneA25;
			tuA25.geneCompartments = cytosol;
			tuA25.compartment = cytosol;

			tuA26.genome = genome;
			tuA26.genes = geneA26;
			tuA26.geneCompartments = cytosol;
			tuA26.compartment = cytosol;

			tuA27.genome = genome;
			tuA27.genes = geneA27;
			tuA27.geneCompartments = cytosol;
			tuA27.compartment = cytosol;

			tuA28.genome = genome;
			tuA28.genes = geneA28;
			tuA28.geneCompartments = cytosol;
			tuA28.compartment = cytosol;

			tuA29.genome = genome;
			tuA29.genes = geneA29;
			tuA29.geneCompartments = cytosol;
			tuA29.compartment = cytosol;

			tuA30.genome = genome;
			tuA30.genes = geneA30;
			tuA30.geneCompartments = cytosol;
			tuA30.compartment = cytosol;

			tuA31.genome = genome;
			tuA31.genes = geneA31;
			tuA31.geneCompartments = cytosol;
			tuA31.compartment = cytosol;

			tuA32.genome = genome;
			tuA32.genes = geneA32;
			tuA32.geneCompartments = cytosol;
			tuA32.compartment = cytosol;

			tuA33.genome = genome;
			tuA33.genes = geneA33;
			tuA33.geneCompartments = cytosol;
			tuA33.compartment = cytosol;

			tuA34.genome = genome;
			tuA34.genes = geneA34;
			tuA34.geneCompartments = cytosol;
			tuA34.compartment = cytosol;

			tuA35.genome = genome;
			tuA35.genes = geneA35;
			tuA35.geneCompartments = cytosol;
			tuA35.compartment = cytosol;

			tuA36.genome = genome;
			tuA36.genes = geneA36;
			tuA36.geneCompartments = cytosol;
			tuA36.compartment = cytosol;

			tuA37.genome = genome;
			tuA37.genes = geneA37;
			tuA37.geneCompartments = cytosol;
			tuA37.compartment = cytosol;

			tuA38.genome = genome;
			tuA38.genes = geneA38;
			tuA38.geneCompartments = cytosol;
			tuA38.compartment = cytosol;

			tuA39.genome = genome;
			tuA39.genes = geneA39;
			tuA39.geneCompartments = cytosol;
			tuA39.compartment = cytosol;

			tuA40.genome = genome;
			tuA40.genes = geneA40;
			tuA40.geneCompartments = cytosol;
			tuA40.compartment = cytosol;

			tuA41.genome = genome;
			tuA41.genes = geneA41;
			tuA41.geneCompartments = cytosol;
			tuA41.compartment = cytosol;

			tuA42.genome = genome;
			tuA42.genes = geneA42;
			tuA42.geneCompartments = cytosol;
			tuA42.compartment = cytosol;

			tuA43.genome = genome;
			tuA43.genes = geneA43;
			tuA43.geneCompartments = cytosol;
			tuA43.compartment = cytosol;

			tuA44.genome = genome;
			tuA44.genes = geneA44;
			tuA44.geneCompartments = cytosol;
			tuA44.compartment = cytosol;

			tuA45.genome = genome;
			tuA45.genes = geneA45;
			tuA45.geneCompartments = cytosol;
			tuA45.compartment = cytosol;

			tuA46.genome = genome;
			tuA46.genes = geneA46;
			tuA46.geneCompartments = cytosol;
			tuA46.compartment = cytosol;

			tuA47.genome = genome;
			tuA47.genes = geneA47;
			tuA47.geneCompartments = cytosol;
			tuA47.compartment = cytosol;

			tuA48.genome = genome;
			tuA48.genes = geneA48;
			tuA48.geneCompartments = cytosol;
			tuA48.compartment = cytosol;

			tuA49.genome = genome;
			tuA49.genes = geneA49;
			tuA49.geneCompartments = cytosol;
			tuA49.compartment = cytosol;

			tuA50.genome = genome;
			tuA50.genes = geneA50;
			tuA50.geneCompartments = cytosol;
			tuA50.compartment = cytosol;

			monA1.gene = geneA1;
			monA1.geneCompartments = cytosol;
			monA1.compartment = cytosol;

			monA2.gene = geneA2;
			monA2.geneCompartments = cytosol;
			monA2.compartment = cytosol;

			monA3.gene = geneA3;
			monA3.geneCompartments = cytosol;
			monA3.compartment = cytosol;

			monA4.gene = geneA4;
			monA4.geneCompartments = cytosol;
			monA4.compartment = cytosol;

			monA5.gene = geneA5;
			monA5.geneCompartments = cytosol;
			monA5.compartment = cytosol;

			monA6.gene = geneA6;
			monA6.geneCompartments = cytosol;
			monA6.compartment = cytosol;

			monA7.gene = geneA7;
			monA7.geneCompartments = cytosol;
			monA7.compartment = cytosol;

			monA8.gene = geneA8;
			monA8.geneCompartments = cytosol;
			monA8.compartment = cytosol;

			monA9.gene = geneA9;
			monA9.geneCompartments = cytosol;
			monA9.compartment = cytosol;

			monA10.gene = geneA10;
			monA10.geneCompartments = cytosol;
			monA10.compartment = cytosol;

			monA11.gene = geneA11;
			monA11.geneCompartments = cytosol;
			monA11.compartment = cytosol;

			monA12.gene = geneA12;
			monA12.geneCompartments = cytosol;
			monA12.compartment = cytosol;

			monA13.gene = geneA13;
			monA13.geneCompartments = cytosol;
			monA13.compartment = cytosol;

			monA14.gene = geneA14;
			monA14.geneCompartments = cytosol;
			monA14.compartment = cytosol;

			monA15.gene = geneA15;
			monA15.geneCompartments = cytosol;
			monA15.compartment = cytosol;

			monA16.gene = geneA16;
			monA16.geneCompartments = cytosol;
			monA16.compartment = cytosol;

			monA17.gene = geneA17;
			monA17.geneCompartments = cytosol;
			monA17.compartment = cytosol;

			monA18.gene = geneA18;
			monA18.geneCompartments = cytosol;
			monA18.compartment = cytosol;

			monA19.gene = geneA19;
			monA19.geneCompartments = cytosol;
			monA19.compartment = cytosol;

			monA20.gene = geneA20;
			monA20.geneCompartments = cytosol;
			monA20.compartment = cytosol;

			monA21.gene = geneA21;
			monA21.geneCompartments = cytosol;
			monA21.compartment = cytosol;

			monA22.gene = geneA22;
			monA22.geneCompartments = cytosol;
			monA22.compartment = cytosol;

			monA23.gene = geneA23;
			monA23.geneCompartments = cytosol;
			monA23.compartment = cytosol;

			monA24.gene = geneA24;
			monA24.geneCompartments = cytosol;
			monA24.compartment = cytosol;

			monA25.gene = geneA25;
			monA25.geneCompartments = cytosol;
			monA25.compartment = cytosol;

			monA26.gene = geneA26;
			monA26.geneCompartments = cytosol;
			monA26.compartment = cytosol;

			monA27.gene = geneA27;
			monA27.geneCompartments = cytosol;
			monA27.compartment = cytosol;

			monA28.gene = geneA28;
			monA28.geneCompartments = cytosol;
			monA28.compartment = cytosol;

			monA29.gene = geneA29;
			monA29.geneCompartments = cytosol;
			monA29.compartment = cytosol;

			monA30.gene = geneA30;
			monA30.geneCompartments = cytosol;
			monA30.compartment = cytosol;

			monA31.gene = geneA31;
			monA31.geneCompartments = cytosol;
			monA31.compartment = cytosol;

			monA32.gene = geneA32;
			monA32.geneCompartments = cytosol;
			monA32.compartment = cytosol;

			monA33.gene = geneA33;
			monA33.geneCompartments = cytosol;
			monA33.compartment = cytosol;

			monA34.gene = geneA34;
			monA34.geneCompartments = cytosol;
			monA34.compartment = cytosol;

			monA35.gene = geneA35;
			monA35.geneCompartments = cytosol;
			monA35.compartment = cytosol;

			monA36.gene = geneA36;
			monA36.geneCompartments = cytosol;
			monA36.compartment = cytosol;

			monA37.gene = geneA37;
			monA37.geneCompartments = cytosol;
			monA37.compartment = cytosol;

			monA38.gene = geneA38;
			monA38.geneCompartments = cytosol;
			monA38.compartment = cytosol;

			monA39.gene = geneA39;
			monA39.geneCompartments = cytosol;
			monA39.compartment = cytosol;

			monA40.gene = geneA40;
			monA40.geneCompartments = cytosol;
			monA40.compartment = cytosol;

			monA41.gene = geneA41;
			monA41.geneCompartments = cytosol;
			monA41.compartment = cytosol;

			monA42.gene = geneA42;
			monA42.geneCompartments = cytosol;
			monA42.compartment = cytosol;

			monA43.gene = geneA43;
			monA43.geneCompartments = cytosol;
			monA43.compartment = cytosol;

			monA44.gene = geneA44;
			monA44.geneCompartments = cytosol;
			monA44.compartment = cytosol;

			monA45.gene = geneA45;
			monA45.geneCompartments = cytosol;
			monA45.compartment = cytosol;

			monA46.gene = geneA46;
			monA46.geneCompartments = cytosol;
			monA46.compartment = cytosol;

			monA47.gene = geneA47;
			monA47.geneCompartments = cytosol;
			monA47.compartment = cytosol;

			monA48.gene = geneA48;
			monA48.geneCompartments = cytosol;
			monA48.compartment = cytosol;

			monA49.gene = geneA49;
			monA49.geneCompartments = cytosol;
			monA49.compartment = cytosol;

			monA50.gene = geneA50;
			monA50.geneCompartments = cytosol;
			monA50.compartment = cytosol;

			cpxAA1.proteinMonomers = monA1;
			cpxAA1.proteinMonomerCompartments = cytosol;
			cpxAA1.proteinMonomerCoefficients = 4;
			cpxAA1.compartment = cytosol;
			cpxAA1.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA2.proteinMonomers = monA2;
			cpxAA2.proteinMonomerCompartments = cytosol;
			cpxAA2.proteinMonomerCoefficients = 4;
			cpxAA2.compartment = cytosol;
			cpxAA2.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA3.proteinMonomers = monA3;
			cpxAA3.proteinMonomerCompartments = cytosol;
			cpxAA3.proteinMonomerCoefficients = 4;
			cpxAA3.compartment = cytosol;
			cpxAA3.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA4.proteinMonomers = monA4;
			cpxAA4.proteinMonomerCompartments = cytosol;
			cpxAA4.proteinMonomerCoefficients = 4;
			cpxAA4.compartment = cytosol;
			cpxAA4.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA5.proteinMonomers = monA5;
			cpxAA5.proteinMonomerCompartments = cytosol;
			cpxAA5.proteinMonomerCoefficients = 4;
			cpxAA5.compartment = cytosol;
			cpxAA5.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA6.proteinMonomers = monA6;
			cpxAA6.proteinMonomerCompartments = cytosol;
			cpxAA6.proteinMonomerCoefficients = 4;
			cpxAA6.compartment = cytosol;
			cpxAA6.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA7.proteinMonomers = monA7;
			cpxAA7.proteinMonomerCompartments = cytosol;
			cpxAA7.proteinMonomerCoefficients = 4;
			cpxAA7.compartment = cytosol;
			cpxAA7.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA8.proteinMonomers = monA8;
			cpxAA8.proteinMonomerCompartments = cytosol;
			cpxAA8.proteinMonomerCoefficients = 4;
			cpxAA8.compartment = cytosol;
			cpxAA8.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA9.proteinMonomers = monA9;
			cpxAA9.proteinMonomerCompartments = cytosol;
			cpxAA9.proteinMonomerCoefficients = 4;
			cpxAA9.compartment = cytosol;
			cpxAA9.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA10.proteinMonomers = monA10;
			cpxAA10.proteinMonomerCompartments = cytosol;
			cpxAA10.proteinMonomerCoefficients = 4;
			cpxAA10.compartment = cytosol;
			cpxAA10.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA11.proteinMonomers = monA11;
			cpxAA11.proteinMonomerCompartments = cytosol;
			cpxAA11.proteinMonomerCoefficients = 4;
			cpxAA11.compartment = cytosol;
			cpxAA11.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA12.proteinMonomers = monA12;
			cpxAA12.proteinMonomerCompartments = cytosol;
			cpxAA12.proteinMonomerCoefficients = 4;
			cpxAA12.compartment = cytosol;
			cpxAA12.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA13.proteinMonomers = monA13;
			cpxAA13.proteinMonomerCompartments = cytosol;
			cpxAA13.proteinMonomerCoefficients = 4;
			cpxAA13.compartment = cytosol;
			cpxAA13.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA14.proteinMonomers = monA14;
			cpxAA14.proteinMonomerCompartments = cytosol;
			cpxAA14.proteinMonomerCoefficients = 4;
			cpxAA14.compartment = cytosol;
			cpxAA14.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA15.proteinMonomers = monA15;
			cpxAA15.proteinMonomerCompartments = cytosol;
			cpxAA15.proteinMonomerCoefficients = 4;
			cpxAA15.compartment = cytosol;
			cpxAA15.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA16.proteinMonomers = monA16;
			cpxAA16.proteinMonomerCompartments = cytosol;
			cpxAA16.proteinMonomerCoefficients = 4;
			cpxAA16.compartment = cytosol;
			cpxAA16.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA17.proteinMonomers = monA17;
			cpxAA17.proteinMonomerCompartments = cytosol;
			cpxAA17.proteinMonomerCoefficients = 4;
			cpxAA17.compartment = cytosol;
			cpxAA17.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA18.proteinMonomers = monA18;
			cpxAA18.proteinMonomerCompartments = cytosol;
			cpxAA18.proteinMonomerCoefficients = 4;
			cpxAA18.compartment = cytosol;
			cpxAA18.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA19.proteinMonomers = monA19;
			cpxAA19.proteinMonomerCompartments = cytosol;
			cpxAA19.proteinMonomerCoefficients = 4;
			cpxAA19.compartment = cytosol;
			cpxAA19.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA20.proteinMonomers = monA20;
			cpxAA20.proteinMonomerCompartments = cytosol;
			cpxAA20.proteinMonomerCoefficients = 4;
			cpxAA20.compartment = cytosol;
			cpxAA20.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA21.proteinMonomers = monA21;
			cpxAA21.proteinMonomerCompartments = cytosol;
			cpxAA21.proteinMonomerCoefficients = 4;
			cpxAA21.compartment = cytosol;
			cpxAA21.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA22.proteinMonomers = monA22;
			cpxAA22.proteinMonomerCompartments = cytosol;
			cpxAA22.proteinMonomerCoefficients = 4;
			cpxAA22.compartment = cytosol;
			cpxAA22.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA23.proteinMonomers = monA23;
			cpxAA23.proteinMonomerCompartments = cytosol;
			cpxAA23.proteinMonomerCoefficients = 4;
			cpxAA23.compartment = cytosol;
			cpxAA23.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA24.proteinMonomers = monA24;
			cpxAA24.proteinMonomerCompartments = cytosol;
			cpxAA24.proteinMonomerCoefficients = 4;
			cpxAA24.compartment = cytosol;
			cpxAA24.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA25.proteinMonomers = monA25;
			cpxAA25.proteinMonomerCompartments = cytosol;
			cpxAA25.proteinMonomerCoefficients = 4;
			cpxAA25.compartment = cytosol;
			cpxAA25.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA26.proteinMonomers = monA26;
			cpxAA26.proteinMonomerCompartments = cytosol;
			cpxAA26.proteinMonomerCoefficients = 4;
			cpxAA26.compartment = cytosol;
			cpxAA26.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA27.proteinMonomers = monA27;
			cpxAA27.proteinMonomerCompartments = cytosol;
			cpxAA27.proteinMonomerCoefficients = 4;
			cpxAA27.compartment = cytosol;
			cpxAA27.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA28.proteinMonomers = monA28;
			cpxAA28.proteinMonomerCompartments = cytosol;
			cpxAA28.proteinMonomerCoefficients = 4;
			cpxAA28.compartment = cytosol;
			cpxAA28.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA29.proteinMonomers = monA29;
			cpxAA29.proteinMonomerCompartments = cytosol;
			cpxAA29.proteinMonomerCoefficients = 4;
			cpxAA29.compartment = cytosol;
			cpxAA29.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA30.proteinMonomers = monA30;
			cpxAA30.proteinMonomerCompartments = cytosol;
			cpxAA30.proteinMonomerCoefficients = 4;
			cpxAA30.compartment = cytosol;
			cpxAA30.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA31.proteinMonomers = monA31;
			cpxAA31.proteinMonomerCompartments = cytosol;
			cpxAA31.proteinMonomerCoefficients = 4;
			cpxAA31.compartment = cytosol;
			cpxAA31.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA32.proteinMonomers = monA32;
			cpxAA32.proteinMonomerCompartments = cytosol;
			cpxAA32.proteinMonomerCoefficients = 4;
			cpxAA32.compartment = cytosol;
			cpxAA32.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA33.proteinMonomers = monA33;
			cpxAA33.proteinMonomerCompartments = cytosol;
			cpxAA33.proteinMonomerCoefficients = 4;
			cpxAA33.compartment = cytosol;
			cpxAA33.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA34.proteinMonomers = monA34;
			cpxAA34.proteinMonomerCompartments = cytosol;
			cpxAA34.proteinMonomerCoefficients = 4;
			cpxAA34.compartment = cytosol;
			cpxAA34.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA35.proteinMonomers = monA35;
			cpxAA35.proteinMonomerCompartments = cytosol;
			cpxAA35.proteinMonomerCoefficients = 4;
			cpxAA35.compartment = cytosol;
			cpxAA35.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA36.proteinMonomers = monA36;
			cpxAA36.proteinMonomerCompartments = cytosol;
			cpxAA36.proteinMonomerCoefficients = 4;
			cpxAA36.compartment = cytosol;
			cpxAA36.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA37.proteinMonomers = monA37;
			cpxAA37.proteinMonomerCompartments = cytosol;
			cpxAA37.proteinMonomerCoefficients = 4;
			cpxAA37.compartment = cytosol;
			cpxAA37.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA38.proteinMonomers = monA38;
			cpxAA38.proteinMonomerCompartments = cytosol;
			cpxAA38.proteinMonomerCoefficients = 4;
			cpxAA38.compartment = cytosol;
			cpxAA38.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA39.proteinMonomers = monA39;
			cpxAA39.proteinMonomerCompartments = cytosol;
			cpxAA39.proteinMonomerCoefficients = 4;
			cpxAA39.compartment = cytosol;
			cpxAA39.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA40.proteinMonomers = monA40;
			cpxAA40.proteinMonomerCompartments = cytosol;
			cpxAA40.proteinMonomerCoefficients = 4;
			cpxAA40.compartment = cytosol;
			cpxAA40.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA41.proteinMonomers = monA41;
			cpxAA41.proteinMonomerCompartments = cytosol;
			cpxAA41.proteinMonomerCoefficients = 4;
			cpxAA41.compartment = cytosol;
			cpxAA41.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA42.proteinMonomers = monA42;
			cpxAA42.proteinMonomerCompartments = cytosol;
			cpxAA42.proteinMonomerCoefficients = 4;
			cpxAA42.compartment = cytosol;
			cpxAA42.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA43.proteinMonomers = monA43;
			cpxAA43.proteinMonomerCompartments = cytosol;
			cpxAA43.proteinMonomerCoefficients = 4;
			cpxAA43.compartment = cytosol;
			cpxAA43.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA44.proteinMonomers = monA44;
			cpxAA44.proteinMonomerCompartments = cytosol;
			cpxAA44.proteinMonomerCoefficients = 4;
			cpxAA44.compartment = cytosol;
			cpxAA44.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA45.proteinMonomers = monA45;
			cpxAA45.proteinMonomerCompartments = cytosol;
			cpxAA45.proteinMonomerCoefficients = 4;
			cpxAA45.compartment = cytosol;
			cpxAA45.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA46.proteinMonomers = monA46;
			cpxAA46.proteinMonomerCompartments = cytosol;
			cpxAA46.proteinMonomerCoefficients = 4;
			cpxAA46.compartment = cytosol;
			cpxAA46.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA47.proteinMonomers = monA47;
			cpxAA47.proteinMonomerCompartments = cytosol;
			cpxAA47.proteinMonomerCoefficients = 4;
			cpxAA47.compartment = cytosol;
			cpxAA47.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA48.proteinMonomers = monA48;
			cpxAA48.proteinMonomerCompartments = cytosol;
			cpxAA48.proteinMonomerCoefficients = 4;
			cpxAA48.compartment = cytosol;
			cpxAA48.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA49.proteinMonomers = monA49;
			cpxAA49.proteinMonomerCompartments = cytosol;
			cpxAA49.proteinMonomerCoefficients = 4;
			cpxAA49.compartment = cytosol;
			cpxAA49.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA50.proteinMonomers = monA50;
			cpxAA50.proteinMonomerCompartments = cytosol;
			cpxAA50.proteinMonomerCoefficients = 4;
			cpxAA50.compartment = cytosol;
			cpxAA50.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			genome.genes = [genome.genes; geneA1; geneA2; geneA3; geneA4; geneA5; geneA6; geneA7; geneA8; geneA9; geneA10; geneA11; geneA12; geneA13; geneA14; geneA15; geneA16; geneA17; geneA18; geneA19; geneA20; geneA21; geneA22; geneA23; geneA24; geneA25; geneA26; geneA27; geneA28; geneA29; geneA30; geneA31; geneA32; geneA33; geneA34; geneA35; geneA36; geneA37; geneA38; geneA39; geneA40; geneA41; geneA42; geneA43; geneA44; geneA45; geneA46; geneA47; geneA48; geneA49; geneA50];
			genome.transcriptionUnits = [genome.transcriptionUnits; tuA1; tuA2; tuA3; tuA4; tuA5; tuA6; tuA7; tuA8; tuA9; tuA10; tuA11; tuA12; tuA13; tuA14; tuA15; tuA16; tuA17; tuA18; tuA19; tuA20; tuA21; tuA22; tuA23; tuA24; tuA25; tuA26; tuA27; tuA28; tuA29; tuA30; tuA31; tuA32; tuA33; tuA34; tuA35; tuA36; tuA37; tuA38; tuA39; tuA40; tuA41; tuA42; tuA43; tuA44; tuA45; tuA46; tuA47; tuA48; tuA49; tuA50];
			kb.genes = genome.genes;
			kb.transcriptionUnits = genome.transcriptionUnits;
			kb.proteinMonomers = [kb.proteinMonomers; monA1; monA2; monA3; monA4; monA5; monA6; monA7; monA8; monA9; monA10; monA11; monA12; monA13; monA14; monA15; monA16; monA17; monA18; monA19; monA20; monA21; monA22; monA23; monA24; monA25; monA26; monA27; monA28; monA29; monA30; monA31; monA32; monA33; monA34; monA35; monA36; monA37; monA38; monA39; monA40; monA41; monA42; monA43; monA44; monA45; monA46; monA47; monA48; monA49; monA50];
			kb.proteinComplexs = [kb.proteinComplexs; cpxAA1; cpxAA2; cpxAA3; cpxAA4; cpxAA5; cpxAA6; cpxAA7; cpxAA8; cpxAA9; cpxAA10; cpxAA11; cpxAA12; cpxAA13; cpxAA14; cpxAA15; cpxAA16; cpxAA17; cpxAA18; cpxAA19; cpxAA20; cpxAA21; cpxAA22; cpxAA23; cpxAA24; cpxAA25; cpxAA26; cpxAA27; cpxAA28; cpxAA29; cpxAA30; cpxAA31; cpxAA32; cpxAA33; cpxAA34; cpxAA35; cpxAA36; cpxAA37; cpxAA38; cpxAA39; cpxAA40; cpxAA41; cpxAA42; cpxAA43; cpxAA44; cpxAA45; cpxAA46; cpxAA47; cpxAA48; cpxAA49; cpxAA50];
			this.modifyNetworkStructure@edu.stanford.covert.cell.sim.runners.SimulationRunner(kb);

		end		function modifyNetworkParameters(~, sim)
			g = sim.gene;
			time = sim.state('Time');
			rna = sim.state('Rna');
			trn = sim.process('Transcription');
			nascentMRNAIndexs = find(any(rna.nascentRNAGeneComposition(g.mRNAIndexs, :), 1));
			[~, modTuIndexs] = ismember({'TuA1', 'TuA2', 'TuA3', 'TuA4', 'TuA5', 'TuA6', 'TuA7', 'TuA8', 'TuA9', 'TuA10', 'TuA11', 'TuA12', 'TuA13', 'TuA14', 'TuA15', 'TuA16', 'TuA17', 'TuA18', 'TuA19', 'TuA20', 'TuA21', 'TuA22', 'TuA23', 'TuA24', 'TuA25', 'TuA26', 'TuA27', 'TuA28', 'TuA29', 'TuA30', 'TuA31', 'TuA32', 'TuA33', 'TuA34', 'TuA35', 'TuA36', 'TuA37', 'TuA38', 'TuA39', 'TuA40', 'TuA41', 'TuA42', 'TuA43', 'TuA44', 'TuA45', 'TuA46', 'TuA47', 'TuA48', 'TuA49', 'TuA50'}, rna.wholeCellModelIDs(rna.nascentIndexs));
			tuBindProb = trn.transcriptionUnitBindingProbabilities;
			meanBindProb = mean(tuBindProb(setdiff(nascentMRNAIndexs, modTuIndexs)));
			rnaDecayRates = rna.decayRates(rna.matureIndexs);
			modTuBindProb = tuBindProb;
			modTuBindProb(modTuIndexs) = meanBindProb;
			modTuBindProb(nascentMRNAIndexs) = modTuBindProb(nascentMRNAIndexs) * sum(tuBindProb(nascentMRNAIndexs)) / sum(modTuBindProb(nascentMRNAIndexs));
			modRnaExp = (rna.nascentRNAMatureRNAComposition * modTuBindProb) ./ (log(2) / time.cellCycleLength + rnaDecayRates);
			modRnaExp = modRnaExp / sum(modRnaExp);
			trn.transcriptionUnitBindingProbabilities = modTuBindProb;
			rna.expression(rna.matureIndexs) = modRnaExp;
		end
	end
end