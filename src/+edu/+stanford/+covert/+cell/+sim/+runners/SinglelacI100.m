classdef SingleGene100 < edu.stanford.covert.cell.sim.runners.SimulationRunner
	methods
		function this = SingleGene100(varargin)
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
			geneA1 = Gene(kb, NaN, {'GeneA1'}, {'Gene A1'}, {'genA1'}, {''}, {'mRNA'}, 0, {''}, 580177, 1083, true, {''}, meanHL, expression);
			geneA2 = Gene(kb, NaN, {'GeneA2'}, {'Gene A2'}, {'genA2'}, {''}, {'mRNA'}, 0, {''}, 581360, 1083, true, {''}, meanHL, expression);
			geneA3 = Gene(kb, NaN, {'GeneA3'}, {'Gene A3'}, {'genA3'}, {''}, {'mRNA'}, 0, {''}, 582543, 1083, true, {''}, meanHL, expression);
			geneA4 = Gene(kb, NaN, {'GeneA4'}, {'Gene A4'}, {'genA4'}, {''}, {'mRNA'}, 0, {''}, 583726, 1083, true, {''}, meanHL, expression);
			geneA5 = Gene(kb, NaN, {'GeneA5'}, {'Gene A5'}, {'genA5'}, {''}, {'mRNA'}, 0, {''}, 584909, 1083, true, {''}, meanHL, expression);
			geneA6 = Gene(kb, NaN, {'GeneA6'}, {'Gene A6'}, {'genA6'}, {''}, {'mRNA'}, 0, {''}, 586092, 1083, true, {''}, meanHL, expression);
			geneA7 = Gene(kb, NaN, {'GeneA7'}, {'Gene A7'}, {'genA7'}, {''}, {'mRNA'}, 0, {''}, 587275, 1083, true, {''}, meanHL, expression);
			geneA8 = Gene(kb, NaN, {'GeneA8'}, {'Gene A8'}, {'genA8'}, {''}, {'mRNA'}, 0, {''}, 588458, 1083, true, {''}, meanHL, expression);
			geneA9 = Gene(kb, NaN, {'GeneA9'}, {'Gene A9'}, {'genA9'}, {''}, {'mRNA'}, 0, {''}, 589641, 1083, true, {''}, meanHL, expression);
			geneA10 = Gene(kb, NaN, {'GeneA10'}, {'Gene A10'}, {'genA10'}, {''}, {'mRNA'}, 0, {''}, 590824, 1083, true, {''}, meanHL, expression);
			geneA11 = Gene(kb, NaN, {'GeneA11'}, {'Gene A11'}, {'genA11'}, {''}, {'mRNA'}, 0, {''}, 592007, 1083, true, {''}, meanHL, expression);
			geneA12 = Gene(kb, NaN, {'GeneA12'}, {'Gene A12'}, {'genA12'}, {''}, {'mRNA'}, 0, {''}, 593190, 1083, true, {''}, meanHL, expression);
			geneA13 = Gene(kb, NaN, {'GeneA13'}, {'Gene A13'}, {'genA13'}, {''}, {'mRNA'}, 0, {''}, 594373, 1083, true, {''}, meanHL, expression);
			geneA14 = Gene(kb, NaN, {'GeneA14'}, {'Gene A14'}, {'genA14'}, {''}, {'mRNA'}, 0, {''}, 595556, 1083, true, {''}, meanHL, expression);
			geneA15 = Gene(kb, NaN, {'GeneA15'}, {'Gene A15'}, {'genA15'}, {''}, {'mRNA'}, 0, {''}, 596739, 1083, true, {''}, meanHL, expression);
			geneA16 = Gene(kb, NaN, {'GeneA16'}, {'Gene A16'}, {'genA16'}, {''}, {'mRNA'}, 0, {''}, 597922, 1083, true, {''}, meanHL, expression);
			geneA17 = Gene(kb, NaN, {'GeneA17'}, {'Gene A17'}, {'genA17'}, {''}, {'mRNA'}, 0, {''}, 599105, 1083, true, {''}, meanHL, expression);
			geneA18 = Gene(kb, NaN, {'GeneA18'}, {'Gene A18'}, {'genA18'}, {''}, {'mRNA'}, 0, {''}, 600288, 1083, true, {''}, meanHL, expression);
			geneA19 = Gene(kb, NaN, {'GeneA19'}, {'Gene A19'}, {'genA19'}, {''}, {'mRNA'}, 0, {''}, 601471, 1083, true, {''}, meanHL, expression);
			geneA20 = Gene(kb, NaN, {'GeneA20'}, {'Gene A20'}, {'genA20'}, {''}, {'mRNA'}, 0, {''}, 602654, 1083, true, {''}, meanHL, expression);
			geneA21 = Gene(kb, NaN, {'GeneA21'}, {'Gene A21'}, {'genA21'}, {''}, {'mRNA'}, 0, {''}, 603837, 1083, true, {''}, meanHL, expression);
			geneA22 = Gene(kb, NaN, {'GeneA22'}, {'Gene A22'}, {'genA22'}, {''}, {'mRNA'}, 0, {''}, 605020, 1083, true, {''}, meanHL, expression);
			geneA23 = Gene(kb, NaN, {'GeneA23'}, {'Gene A23'}, {'genA23'}, {''}, {'mRNA'}, 0, {''}, 606203, 1083, true, {''}, meanHL, expression);
			geneA24 = Gene(kb, NaN, {'GeneA24'}, {'Gene A24'}, {'genA24'}, {''}, {'mRNA'}, 0, {''}, 607386, 1083, true, {''}, meanHL, expression);
			geneA25 = Gene(kb, NaN, {'GeneA25'}, {'Gene A25'}, {'genA25'}, {''}, {'mRNA'}, 0, {''}, 608569, 1083, true, {''}, meanHL, expression);
			geneA26 = Gene(kb, NaN, {'GeneA26'}, {'Gene A26'}, {'genA26'}, {''}, {'mRNA'}, 0, {''}, 609752, 1083, true, {''}, meanHL, expression);
			geneA27 = Gene(kb, NaN, {'GeneA27'}, {'Gene A27'}, {'genA27'}, {''}, {'mRNA'}, 0, {''}, 610935, 1083, true, {''}, meanHL, expression);
			geneA28 = Gene(kb, NaN, {'GeneA28'}, {'Gene A28'}, {'genA28'}, {''}, {'mRNA'}, 0, {''}, 612118, 1083, true, {''}, meanHL, expression);
			geneA29 = Gene(kb, NaN, {'GeneA29'}, {'Gene A29'}, {'genA29'}, {''}, {'mRNA'}, 0, {''}, 613301, 1083, true, {''}, meanHL, expression);
			geneA30 = Gene(kb, NaN, {'GeneA30'}, {'Gene A30'}, {'genA30'}, {''}, {'mRNA'}, 0, {''}, 614484, 1083, true, {''}, meanHL, expression);
			geneA31 = Gene(kb, NaN, {'GeneA31'}, {'Gene A31'}, {'genA31'}, {''}, {'mRNA'}, 0, {''}, 615667, 1083, true, {''}, meanHL, expression);
			geneA32 = Gene(kb, NaN, {'GeneA32'}, {'Gene A32'}, {'genA32'}, {''}, {'mRNA'}, 0, {''}, 616850, 1083, true, {''}, meanHL, expression);
			geneA33 = Gene(kb, NaN, {'GeneA33'}, {'Gene A33'}, {'genA33'}, {''}, {'mRNA'}, 0, {''}, 618033, 1083, true, {''}, meanHL, expression);
			geneA34 = Gene(kb, NaN, {'GeneA34'}, {'Gene A34'}, {'genA34'}, {''}, {'mRNA'}, 0, {''}, 619216, 1083, true, {''}, meanHL, expression);
			geneA35 = Gene(kb, NaN, {'GeneA35'}, {'Gene A35'}, {'genA35'}, {''}, {'mRNA'}, 0, {''}, 620399, 1083, true, {''}, meanHL, expression);
			geneA36 = Gene(kb, NaN, {'GeneA36'}, {'Gene A36'}, {'genA36'}, {''}, {'mRNA'}, 0, {''}, 621582, 1083, true, {''}, meanHL, expression);
			geneA37 = Gene(kb, NaN, {'GeneA37'}, {'Gene A37'}, {'genA37'}, {''}, {'mRNA'}, 0, {''}, 622765, 1083, true, {''}, meanHL, expression);
			geneA38 = Gene(kb, NaN, {'GeneA38'}, {'Gene A38'}, {'genA38'}, {''}, {'mRNA'}, 0, {''}, 623948, 1083, true, {''}, meanHL, expression);
			geneA39 = Gene(kb, NaN, {'GeneA39'}, {'Gene A39'}, {'genA39'}, {''}, {'mRNA'}, 0, {''}, 625131, 1083, true, {''}, meanHL, expression);
			geneA40 = Gene(kb, NaN, {'GeneA40'}, {'Gene A40'}, {'genA40'}, {''}, {'mRNA'}, 0, {''}, 626314, 1083, true, {''}, meanHL, expression);
			geneA41 = Gene(kb, NaN, {'GeneA41'}, {'Gene A41'}, {'genA41'}, {''}, {'mRNA'}, 0, {''}, 627497, 1083, true, {''}, meanHL, expression);
			geneA42 = Gene(kb, NaN, {'GeneA42'}, {'Gene A42'}, {'genA42'}, {''}, {'mRNA'}, 0, {''}, 628680, 1083, true, {''}, meanHL, expression);
			geneA43 = Gene(kb, NaN, {'GeneA43'}, {'Gene A43'}, {'genA43'}, {''}, {'mRNA'}, 0, {''}, 629863, 1083, true, {''}, meanHL, expression);
			geneA44 = Gene(kb, NaN, {'GeneA44'}, {'Gene A44'}, {'genA44'}, {''}, {'mRNA'}, 0, {''}, 631046, 1083, true, {''}, meanHL, expression);
			geneA45 = Gene(kb, NaN, {'GeneA45'}, {'Gene A45'}, {'genA45'}, {''}, {'mRNA'}, 0, {''}, 632229, 1083, true, {''}, meanHL, expression);
			geneA46 = Gene(kb, NaN, {'GeneA46'}, {'Gene A46'}, {'genA46'}, {''}, {'mRNA'}, 0, {''}, 633412, 1083, true, {''}, meanHL, expression);
			geneA47 = Gene(kb, NaN, {'GeneA47'}, {'Gene A47'}, {'genA47'}, {''}, {'mRNA'}, 0, {''}, 634595, 1083, true, {''}, meanHL, expression);
			geneA48 = Gene(kb, NaN, {'GeneA48'}, {'Gene A48'}, {'genA48'}, {''}, {'mRNA'}, 0, {''}, 635778, 1083, true, {''}, meanHL, expression);
			geneA49 = Gene(kb, NaN, {'GeneA49'}, {'Gene A49'}, {'genA49'}, {''}, {'mRNA'}, 0, {''}, 636961, 1083, true, {''}, meanHL, expression);
			geneA50 = Gene(kb, NaN, {'GeneA50'}, {'Gene A50'}, {'genA50'}, {''}, {'mRNA'}, 0, {''}, 638144, 1083, true, {''}, meanHL, expression);
			geneA51 = Gene(kb, NaN, {'GeneA51'}, {'Gene A51'}, {'genA51'}, {''}, {'mRNA'}, 0, {''}, 639327, 1083, true, {''}, meanHL, expression);
			geneA52 = Gene(kb, NaN, {'GeneA52'}, {'Gene A52'}, {'genA52'}, {''}, {'mRNA'}, 0, {''}, 640510, 1083, true, {''}, meanHL, expression);
			geneA53 = Gene(kb, NaN, {'GeneA53'}, {'Gene A53'}, {'genA53'}, {''}, {'mRNA'}, 0, {''}, 641693, 1083, true, {''}, meanHL, expression);
			geneA54 = Gene(kb, NaN, {'GeneA54'}, {'Gene A54'}, {'genA54'}, {''}, {'mRNA'}, 0, {''}, 642876, 1083, true, {''}, meanHL, expression);
			geneA55 = Gene(kb, NaN, {'GeneA55'}, {'Gene A55'}, {'genA55'}, {''}, {'mRNA'}, 0, {''}, 644059, 1083, true, {''}, meanHL, expression);
			geneA56 = Gene(kb, NaN, {'GeneA56'}, {'Gene A56'}, {'genA56'}, {''}, {'mRNA'}, 0, {''}, 645242, 1083, true, {''}, meanHL, expression);
			geneA57 = Gene(kb, NaN, {'GeneA57'}, {'Gene A57'}, {'genA57'}, {''}, {'mRNA'}, 0, {''}, 646425, 1083, true, {''}, meanHL, expression);
			geneA58 = Gene(kb, NaN, {'GeneA58'}, {'Gene A58'}, {'genA58'}, {''}, {'mRNA'}, 0, {''}, 647608, 1083, true, {''}, meanHL, expression);
			geneA59 = Gene(kb, NaN, {'GeneA59'}, {'Gene A59'}, {'genA59'}, {''}, {'mRNA'}, 0, {''}, 648791, 1083, true, {''}, meanHL, expression);
			geneA60 = Gene(kb, NaN, {'GeneA60'}, {'Gene A60'}, {'genA60'}, {''}, {'mRNA'}, 0, {''}, 649974, 1083, true, {''}, meanHL, expression);
			geneA61 = Gene(kb, NaN, {'GeneA61'}, {'Gene A61'}, {'genA61'}, {''}, {'mRNA'}, 0, {''}, 651157, 1083, true, {''}, meanHL, expression);
			geneA62 = Gene(kb, NaN, {'GeneA62'}, {'Gene A62'}, {'genA62'}, {''}, {'mRNA'}, 0, {''}, 652340, 1083, true, {''}, meanHL, expression);
			geneA63 = Gene(kb, NaN, {'GeneA63'}, {'Gene A63'}, {'genA63'}, {''}, {'mRNA'}, 0, {''}, 653523, 1083, true, {''}, meanHL, expression);
			geneA64 = Gene(kb, NaN, {'GeneA64'}, {'Gene A64'}, {'genA64'}, {''}, {'mRNA'}, 0, {''}, 654706, 1083, true, {''}, meanHL, expression);
			geneA65 = Gene(kb, NaN, {'GeneA65'}, {'Gene A65'}, {'genA65'}, {''}, {'mRNA'}, 0, {''}, 655889, 1083, true, {''}, meanHL, expression);
			geneA66 = Gene(kb, NaN, {'GeneA66'}, {'Gene A66'}, {'genA66'}, {''}, {'mRNA'}, 0, {''}, 657072, 1083, true, {''}, meanHL, expression);
			geneA67 = Gene(kb, NaN, {'GeneA67'}, {'Gene A67'}, {'genA67'}, {''}, {'mRNA'}, 0, {''}, 658255, 1083, true, {''}, meanHL, expression);
			geneA68 = Gene(kb, NaN, {'GeneA68'}, {'Gene A68'}, {'genA68'}, {''}, {'mRNA'}, 0, {''}, 659438, 1083, true, {''}, meanHL, expression);
			geneA69 = Gene(kb, NaN, {'GeneA69'}, {'Gene A69'}, {'genA69'}, {''}, {'mRNA'}, 0, {''}, 660621, 1083, true, {''}, meanHL, expression);
			geneA70 = Gene(kb, NaN, {'GeneA70'}, {'Gene A70'}, {'genA70'}, {''}, {'mRNA'}, 0, {''}, 661804, 1083, true, {''}, meanHL, expression);
			geneA71 = Gene(kb, NaN, {'GeneA71'}, {'Gene A71'}, {'genA71'}, {''}, {'mRNA'}, 0, {''}, 662987, 1083, true, {''}, meanHL, expression);
			geneA72 = Gene(kb, NaN, {'GeneA72'}, {'Gene A72'}, {'genA72'}, {''}, {'mRNA'}, 0, {''}, 664170, 1083, true, {''}, meanHL, expression);
			geneA73 = Gene(kb, NaN, {'GeneA73'}, {'Gene A73'}, {'genA73'}, {''}, {'mRNA'}, 0, {''}, 665353, 1083, true, {''}, meanHL, expression);
			geneA74 = Gene(kb, NaN, {'GeneA74'}, {'Gene A74'}, {'genA74'}, {''}, {'mRNA'}, 0, {''}, 666536, 1083, true, {''}, meanHL, expression);
			geneA75 = Gene(kb, NaN, {'GeneA75'}, {'Gene A75'}, {'genA75'}, {''}, {'mRNA'}, 0, {''}, 667719, 1083, true, {''}, meanHL, expression);
			geneA76 = Gene(kb, NaN, {'GeneA76'}, {'Gene A76'}, {'genA76'}, {''}, {'mRNA'}, 0, {''}, 668902, 1083, true, {''}, meanHL, expression);
			geneA77 = Gene(kb, NaN, {'GeneA77'}, {'Gene A77'}, {'genA77'}, {''}, {'mRNA'}, 0, {''}, 670085, 1083, true, {''}, meanHL, expression);
			geneA78 = Gene(kb, NaN, {'GeneA78'}, {'Gene A78'}, {'genA78'}, {''}, {'mRNA'}, 0, {''}, 671268, 1083, true, {''}, meanHL, expression);
			geneA79 = Gene(kb, NaN, {'GeneA79'}, {'Gene A79'}, {'genA79'}, {''}, {'mRNA'}, 0, {''}, 672451, 1083, true, {''}, meanHL, expression);
			geneA80 = Gene(kb, NaN, {'GeneA80'}, {'Gene A80'}, {'genA80'}, {''}, {'mRNA'}, 0, {''}, 673634, 1083, true, {''}, meanHL, expression);
			geneA81 = Gene(kb, NaN, {'GeneA81'}, {'Gene A81'}, {'genA81'}, {''}, {'mRNA'}, 0, {''}, 674817, 1083, true, {''}, meanHL, expression);
			geneA82 = Gene(kb, NaN, {'GeneA82'}, {'Gene A82'}, {'genA82'}, {''}, {'mRNA'}, 0, {''}, 676000, 1083, true, {''}, meanHL, expression);
			geneA83 = Gene(kb, NaN, {'GeneA83'}, {'Gene A83'}, {'genA83'}, {''}, {'mRNA'}, 0, {''}, 677183, 1083, true, {''}, meanHL, expression);
			geneA84 = Gene(kb, NaN, {'GeneA84'}, {'Gene A84'}, {'genA84'}, {''}, {'mRNA'}, 0, {''}, 678366, 1083, true, {''}, meanHL, expression);
			geneA85 = Gene(kb, NaN, {'GeneA85'}, {'Gene A85'}, {'genA85'}, {''}, {'mRNA'}, 0, {''}, 679549, 1083, true, {''}, meanHL, expression);
			geneA86 = Gene(kb, NaN, {'GeneA86'}, {'Gene A86'}, {'genA86'}, {''}, {'mRNA'}, 0, {''}, 680732, 1083, true, {''}, meanHL, expression);
			geneA87 = Gene(kb, NaN, {'GeneA87'}, {'Gene A87'}, {'genA87'}, {''}, {'mRNA'}, 0, {''}, 681915, 1083, true, {''}, meanHL, expression);
			geneA88 = Gene(kb, NaN, {'GeneA88'}, {'Gene A88'}, {'genA88'}, {''}, {'mRNA'}, 0, {''}, 683098, 1083, true, {''}, meanHL, expression);
			geneA89 = Gene(kb, NaN, {'GeneA89'}, {'Gene A89'}, {'genA89'}, {''}, {'mRNA'}, 0, {''}, 684281, 1083, true, {''}, meanHL, expression);
			geneA90 = Gene(kb, NaN, {'GeneA90'}, {'Gene A90'}, {'genA90'}, {''}, {'mRNA'}, 0, {''}, 685464, 1083, true, {''}, meanHL, expression);
			geneA91 = Gene(kb, NaN, {'GeneA91'}, {'Gene A91'}, {'genA91'}, {''}, {'mRNA'}, 0, {''}, 686647, 1083, true, {''}, meanHL, expression);
			geneA92 = Gene(kb, NaN, {'GeneA92'}, {'Gene A92'}, {'genA92'}, {''}, {'mRNA'}, 0, {''}, 687830, 1083, true, {''}, meanHL, expression);
			geneA93 = Gene(kb, NaN, {'GeneA93'}, {'Gene A93'}, {'genA93'}, {''}, {'mRNA'}, 0, {''}, 689013, 1083, true, {''}, meanHL, expression);
			geneA94 = Gene(kb, NaN, {'GeneA94'}, {'Gene A94'}, {'genA94'}, {''}, {'mRNA'}, 0, {''}, 690196, 1083, true, {''}, meanHL, expression);
			geneA95 = Gene(kb, NaN, {'GeneA95'}, {'Gene A95'}, {'genA95'}, {''}, {'mRNA'}, 0, {''}, 691379, 1083, true, {''}, meanHL, expression);
			geneA96 = Gene(kb, NaN, {'GeneA96'}, {'Gene A96'}, {'genA96'}, {''}, {'mRNA'}, 0, {''}, 692562, 1083, true, {''}, meanHL, expression);
			geneA97 = Gene(kb, NaN, {'GeneA97'}, {'Gene A97'}, {'genA97'}, {''}, {'mRNA'}, 0, {''}, 693745, 1083, true, {''}, meanHL, expression);
			geneA98 = Gene(kb, NaN, {'GeneA98'}, {'Gene A98'}, {'genA98'}, {''}, {'mRNA'}, 0, {''}, 694928, 1083, true, {''}, meanHL, expression);
			geneA99 = Gene(kb, NaN, {'GeneA99'}, {'Gene A99'}, {'genA99'}, {''}, {'mRNA'}, 0, {''}, 696111, 1083, true, {''}, meanHL, expression);
			geneA100 = Gene(kb, NaN, {'GeneA100'}, {'Gene A100'}, {'genA100'}, {''}, {'mRNA'}, 0, {''}, 697294, 1083, true, {''}, meanHL, expression);
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
			tuA51 = TranscriptionUnit(kb, NaN, {'TuA51'}, {'Transcription unit A51'}, -35, 6, -35, 6, -1);
			tuA52 = TranscriptionUnit(kb, NaN, {'TuA52'}, {'Transcription unit A52'}, -35, 6, -35, 6, -1);
			tuA53 = TranscriptionUnit(kb, NaN, {'TuA53'}, {'Transcription unit A53'}, -35, 6, -35, 6, -1);
			tuA54 = TranscriptionUnit(kb, NaN, {'TuA54'}, {'Transcription unit A54'}, -35, 6, -35, 6, -1);
			tuA55 = TranscriptionUnit(kb, NaN, {'TuA55'}, {'Transcription unit A55'}, -35, 6, -35, 6, -1);
			tuA56 = TranscriptionUnit(kb, NaN, {'TuA56'}, {'Transcription unit A56'}, -35, 6, -35, 6, -1);
			tuA57 = TranscriptionUnit(kb, NaN, {'TuA57'}, {'Transcription unit A57'}, -35, 6, -35, 6, -1);
			tuA58 = TranscriptionUnit(kb, NaN, {'TuA58'}, {'Transcription unit A58'}, -35, 6, -35, 6, -1);
			tuA59 = TranscriptionUnit(kb, NaN, {'TuA59'}, {'Transcription unit A59'}, -35, 6, -35, 6, -1);
			tuA60 = TranscriptionUnit(kb, NaN, {'TuA60'}, {'Transcription unit A60'}, -35, 6, -35, 6, -1);
			tuA61 = TranscriptionUnit(kb, NaN, {'TuA61'}, {'Transcription unit A61'}, -35, 6, -35, 6, -1);
			tuA62 = TranscriptionUnit(kb, NaN, {'TuA62'}, {'Transcription unit A62'}, -35, 6, -35, 6, -1);
			tuA63 = TranscriptionUnit(kb, NaN, {'TuA63'}, {'Transcription unit A63'}, -35, 6, -35, 6, -1);
			tuA64 = TranscriptionUnit(kb, NaN, {'TuA64'}, {'Transcription unit A64'}, -35, 6, -35, 6, -1);
			tuA65 = TranscriptionUnit(kb, NaN, {'TuA65'}, {'Transcription unit A65'}, -35, 6, -35, 6, -1);
			tuA66 = TranscriptionUnit(kb, NaN, {'TuA66'}, {'Transcription unit A66'}, -35, 6, -35, 6, -1);
			tuA67 = TranscriptionUnit(kb, NaN, {'TuA67'}, {'Transcription unit A67'}, -35, 6, -35, 6, -1);
			tuA68 = TranscriptionUnit(kb, NaN, {'TuA68'}, {'Transcription unit A68'}, -35, 6, -35, 6, -1);
			tuA69 = TranscriptionUnit(kb, NaN, {'TuA69'}, {'Transcription unit A69'}, -35, 6, -35, 6, -1);
			tuA70 = TranscriptionUnit(kb, NaN, {'TuA70'}, {'Transcription unit A70'}, -35, 6, -35, 6, -1);
			tuA71 = TranscriptionUnit(kb, NaN, {'TuA71'}, {'Transcription unit A71'}, -35, 6, -35, 6, -1);
			tuA72 = TranscriptionUnit(kb, NaN, {'TuA72'}, {'Transcription unit A72'}, -35, 6, -35, 6, -1);
			tuA73 = TranscriptionUnit(kb, NaN, {'TuA73'}, {'Transcription unit A73'}, -35, 6, -35, 6, -1);
			tuA74 = TranscriptionUnit(kb, NaN, {'TuA74'}, {'Transcription unit A74'}, -35, 6, -35, 6, -1);
			tuA75 = TranscriptionUnit(kb, NaN, {'TuA75'}, {'Transcription unit A75'}, -35, 6, -35, 6, -1);
			tuA76 = TranscriptionUnit(kb, NaN, {'TuA76'}, {'Transcription unit A76'}, -35, 6, -35, 6, -1);
			tuA77 = TranscriptionUnit(kb, NaN, {'TuA77'}, {'Transcription unit A77'}, -35, 6, -35, 6, -1);
			tuA78 = TranscriptionUnit(kb, NaN, {'TuA78'}, {'Transcription unit A78'}, -35, 6, -35, 6, -1);
			tuA79 = TranscriptionUnit(kb, NaN, {'TuA79'}, {'Transcription unit A79'}, -35, 6, -35, 6, -1);
			tuA80 = TranscriptionUnit(kb, NaN, {'TuA80'}, {'Transcription unit A80'}, -35, 6, -35, 6, -1);
			tuA81 = TranscriptionUnit(kb, NaN, {'TuA81'}, {'Transcription unit A81'}, -35, 6, -35, 6, -1);
			tuA82 = TranscriptionUnit(kb, NaN, {'TuA82'}, {'Transcription unit A82'}, -35, 6, -35, 6, -1);
			tuA83 = TranscriptionUnit(kb, NaN, {'TuA83'}, {'Transcription unit A83'}, -35, 6, -35, 6, -1);
			tuA84 = TranscriptionUnit(kb, NaN, {'TuA84'}, {'Transcription unit A84'}, -35, 6, -35, 6, -1);
			tuA85 = TranscriptionUnit(kb, NaN, {'TuA85'}, {'Transcription unit A85'}, -35, 6, -35, 6, -1);
			tuA86 = TranscriptionUnit(kb, NaN, {'TuA86'}, {'Transcription unit A86'}, -35, 6, -35, 6, -1);
			tuA87 = TranscriptionUnit(kb, NaN, {'TuA87'}, {'Transcription unit A87'}, -35, 6, -35, 6, -1);
			tuA88 = TranscriptionUnit(kb, NaN, {'TuA88'}, {'Transcription unit A88'}, -35, 6, -35, 6, -1);
			tuA89 = TranscriptionUnit(kb, NaN, {'TuA89'}, {'Transcription unit A89'}, -35, 6, -35, 6, -1);
			tuA90 = TranscriptionUnit(kb, NaN, {'TuA90'}, {'Transcription unit A90'}, -35, 6, -35, 6, -1);
			tuA91 = TranscriptionUnit(kb, NaN, {'TuA91'}, {'Transcription unit A91'}, -35, 6, -35, 6, -1);
			tuA92 = TranscriptionUnit(kb, NaN, {'TuA92'}, {'Transcription unit A92'}, -35, 6, -35, 6, -1);
			tuA93 = TranscriptionUnit(kb, NaN, {'TuA93'}, {'Transcription unit A93'}, -35, 6, -35, 6, -1);
			tuA94 = TranscriptionUnit(kb, NaN, {'TuA94'}, {'Transcription unit A94'}, -35, 6, -35, 6, -1);
			tuA95 = TranscriptionUnit(kb, NaN, {'TuA95'}, {'Transcription unit A95'}, -35, 6, -35, 6, -1);
			tuA96 = TranscriptionUnit(kb, NaN, {'TuA96'}, {'Transcription unit A96'}, -35, 6, -35, 6, -1);
			tuA97 = TranscriptionUnit(kb, NaN, {'TuA97'}, {'Transcription unit A97'}, -35, 6, -35, 6, -1);
			tuA98 = TranscriptionUnit(kb, NaN, {'TuA98'}, {'Transcription unit A98'}, -35, 6, -35, 6, -1);
			tuA99 = TranscriptionUnit(kb, NaN, {'TuA99'}, {'Transcription unit A99'}, -35, 6, -35, 6, -1);
			tuA100 = TranscriptionUnit(kb, NaN, {'TuA100'}, {'Transcription unit A100'}, -35, 6, -35, 6, -1);
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
			monA51 = ProteinMonomer(kb, NaN, {'MonomerA51'}, {'Protein monomer A51'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA52 = ProteinMonomer(kb, NaN, {'MonomerA52'}, {'Protein monomer A52'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA53 = ProteinMonomer(kb, NaN, {'MonomerA53'}, {'Protein monomer A53'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA54 = ProteinMonomer(kb, NaN, {'MonomerA54'}, {'Protein monomer A54'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA55 = ProteinMonomer(kb, NaN, {'MonomerA55'}, {'Protein monomer A55'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA56 = ProteinMonomer(kb, NaN, {'MonomerA56'}, {'Protein monomer A56'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA57 = ProteinMonomer(kb, NaN, {'MonomerA57'}, {'Protein monomer A57'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA58 = ProteinMonomer(kb, NaN, {'MonomerA58'}, {'Protein monomer A58'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA59 = ProteinMonomer(kb, NaN, {'MonomerA59'}, {'Protein monomer A59'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA60 = ProteinMonomer(kb, NaN, {'MonomerA60'}, {'Protein monomer A60'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA61 = ProteinMonomer(kb, NaN, {'MonomerA61'}, {'Protein monomer A61'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA62 = ProteinMonomer(kb, NaN, {'MonomerA62'}, {'Protein monomer A62'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA63 = ProteinMonomer(kb, NaN, {'MonomerA63'}, {'Protein monomer A63'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA64 = ProteinMonomer(kb, NaN, {'MonomerA64'}, {'Protein monomer A64'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA65 = ProteinMonomer(kb, NaN, {'MonomerA65'}, {'Protein monomer A65'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA66 = ProteinMonomer(kb, NaN, {'MonomerA66'}, {'Protein monomer A66'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA67 = ProteinMonomer(kb, NaN, {'MonomerA67'}, {'Protein monomer A67'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA68 = ProteinMonomer(kb, NaN, {'MonomerA68'}, {'Protein monomer A68'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA69 = ProteinMonomer(kb, NaN, {'MonomerA69'}, {'Protein monomer A69'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA70 = ProteinMonomer(kb, NaN, {'MonomerA70'}, {'Protein monomer A70'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA71 = ProteinMonomer(kb, NaN, {'MonomerA71'}, {'Protein monomer A71'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA72 = ProteinMonomer(kb, NaN, {'MonomerA72'}, {'Protein monomer A72'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA73 = ProteinMonomer(kb, NaN, {'MonomerA73'}, {'Protein monomer A73'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA74 = ProteinMonomer(kb, NaN, {'MonomerA74'}, {'Protein monomer A74'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA75 = ProteinMonomer(kb, NaN, {'MonomerA75'}, {'Protein monomer A75'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA76 = ProteinMonomer(kb, NaN, {'MonomerA76'}, {'Protein monomer A76'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA77 = ProteinMonomer(kb, NaN, {'MonomerA77'}, {'Protein monomer A77'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA78 = ProteinMonomer(kb, NaN, {'MonomerA78'}, {'Protein monomer A78'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA79 = ProteinMonomer(kb, NaN, {'MonomerA79'}, {'Protein monomer A79'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA80 = ProteinMonomer(kb, NaN, {'MonomerA80'}, {'Protein monomer A80'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA81 = ProteinMonomer(kb, NaN, {'MonomerA81'}, {'Protein monomer A81'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA82 = ProteinMonomer(kb, NaN, {'MonomerA82'}, {'Protein monomer A82'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA83 = ProteinMonomer(kb, NaN, {'MonomerA83'}, {'Protein monomer A83'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA84 = ProteinMonomer(kb, NaN, {'MonomerA84'}, {'Protein monomer A84'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA85 = ProteinMonomer(kb, NaN, {'MonomerA85'}, {'Protein monomer A85'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA86 = ProteinMonomer(kb, NaN, {'MonomerA86'}, {'Protein monomer A86'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA87 = ProteinMonomer(kb, NaN, {'MonomerA87'}, {'Protein monomer A87'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA88 = ProteinMonomer(kb, NaN, {'MonomerA88'}, {'Protein monomer A88'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA89 = ProteinMonomer(kb, NaN, {'MonomerA89'}, {'Protein monomer A89'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA90 = ProteinMonomer(kb, NaN, {'MonomerA90'}, {'Protein monomer A90'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA91 = ProteinMonomer(kb, NaN, {'MonomerA91'}, {'Protein monomer A91'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA92 = ProteinMonomer(kb, NaN, {'MonomerA92'}, {'Protein monomer A92'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA93 = ProteinMonomer(kb, NaN, {'MonomerA93'}, {'Protein monomer A93'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA94 = ProteinMonomer(kb, NaN, {'MonomerA94'}, {'Protein monomer A94'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA95 = ProteinMonomer(kb, NaN, {'MonomerA95'}, {'Protein monomer A95'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA96 = ProteinMonomer(kb, NaN, {'MonomerA96'}, {'Protein monomer A96'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA97 = ProteinMonomer(kb, NaN, {'MonomerA97'}, {'Protein monomer A97'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA98 = ProteinMonomer(kb, NaN, {'MonomerA98'}, {'Protein monomer A98'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA99 = ProteinMonomer(kb, NaN, {'MonomerA99'}, {'Protein monomer A99'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
			monA100 = ProteinMonomer(kb, NaN, {'MonomerA100'}, {'Protein monomer A100'}, {''}, {''}, {''}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, false, {''}, {''}, 0, {''});
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
			cpxAA51 = ProteinComplex(kb, NaN, {'ComplexAA51'}, {'Macromolecular complex AA51'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA52 = ProteinComplex(kb, NaN, {'ComplexAA52'}, {'Macromolecular complex AA52'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA53 = ProteinComplex(kb, NaN, {'ComplexAA53'}, {'Macromolecular complex AA53'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA54 = ProteinComplex(kb, NaN, {'ComplexAA54'}, {'Macromolecular complex AA54'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA55 = ProteinComplex(kb, NaN, {'ComplexAA55'}, {'Macromolecular complex AA55'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA56 = ProteinComplex(kb, NaN, {'ComplexAA56'}, {'Macromolecular complex AA56'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA57 = ProteinComplex(kb, NaN, {'ComplexAA57'}, {'Macromolecular complex AA57'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA58 = ProteinComplex(kb, NaN, {'ComplexAA58'}, {'Macromolecular complex AA58'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA59 = ProteinComplex(kb, NaN, {'ComplexAA59'}, {'Macromolecular complex AA59'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA60 = ProteinComplex(kb, NaN, {'ComplexAA60'}, {'Macromolecular complex AA60'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA61 = ProteinComplex(kb, NaN, {'ComplexAA61'}, {'Macromolecular complex AA61'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA62 = ProteinComplex(kb, NaN, {'ComplexAA62'}, {'Macromolecular complex AA62'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA63 = ProteinComplex(kb, NaN, {'ComplexAA63'}, {'Macromolecular complex AA63'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA64 = ProteinComplex(kb, NaN, {'ComplexAA64'}, {'Macromolecular complex AA64'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA65 = ProteinComplex(kb, NaN, {'ComplexAA65'}, {'Macromolecular complex AA65'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA66 = ProteinComplex(kb, NaN, {'ComplexAA66'}, {'Macromolecular complex AA66'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA67 = ProteinComplex(kb, NaN, {'ComplexAA67'}, {'Macromolecular complex AA67'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA68 = ProteinComplex(kb, NaN, {'ComplexAA68'}, {'Macromolecular complex AA68'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA69 = ProteinComplex(kb, NaN, {'ComplexAA69'}, {'Macromolecular complex AA69'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA70 = ProteinComplex(kb, NaN, {'ComplexAA70'}, {'Macromolecular complex AA70'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA71 = ProteinComplex(kb, NaN, {'ComplexAA71'}, {'Macromolecular complex AA71'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA72 = ProteinComplex(kb, NaN, {'ComplexAA72'}, {'Macromolecular complex AA72'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA73 = ProteinComplex(kb, NaN, {'ComplexAA73'}, {'Macromolecular complex AA73'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA74 = ProteinComplex(kb, NaN, {'ComplexAA74'}, {'Macromolecular complex AA74'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA75 = ProteinComplex(kb, NaN, {'ComplexAA75'}, {'Macromolecular complex AA75'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA76 = ProteinComplex(kb, NaN, {'ComplexAA76'}, {'Macromolecular complex AA76'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA77 = ProteinComplex(kb, NaN, {'ComplexAA77'}, {'Macromolecular complex AA77'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA78 = ProteinComplex(kb, NaN, {'ComplexAA78'}, {'Macromolecular complex AA78'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA79 = ProteinComplex(kb, NaN, {'ComplexAA79'}, {'Macromolecular complex AA79'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA80 = ProteinComplex(kb, NaN, {'ComplexAA80'}, {'Macromolecular complex AA80'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA81 = ProteinComplex(kb, NaN, {'ComplexAA81'}, {'Macromolecular complex AA81'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA82 = ProteinComplex(kb, NaN, {'ComplexAA82'}, {'Macromolecular complex AA82'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA83 = ProteinComplex(kb, NaN, {'ComplexAA83'}, {'Macromolecular complex AA83'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA84 = ProteinComplex(kb, NaN, {'ComplexAA84'}, {'Macromolecular complex AA84'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA85 = ProteinComplex(kb, NaN, {'ComplexAA85'}, {'Macromolecular complex AA85'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA86 = ProteinComplex(kb, NaN, {'ComplexAA86'}, {'Macromolecular complex AA86'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA87 = ProteinComplex(kb, NaN, {'ComplexAA87'}, {'Macromolecular complex AA87'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA88 = ProteinComplex(kb, NaN, {'ComplexAA88'}, {'Macromolecular complex AA88'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA89 = ProteinComplex(kb, NaN, {'ComplexAA89'}, {'Macromolecular complex AA89'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA90 = ProteinComplex(kb, NaN, {'ComplexAA90'}, {'Macromolecular complex AA90'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA91 = ProteinComplex(kb, NaN, {'ComplexAA91'}, {'Macromolecular complex AA91'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA92 = ProteinComplex(kb, NaN, {'ComplexAA92'}, {'Macromolecular complex AA92'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA93 = ProteinComplex(kb, NaN, {'ComplexAA93'}, {'Macromolecular complex AA93'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA94 = ProteinComplex(kb, NaN, {'ComplexAA94'}, {'Macromolecular complex AA94'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA95 = ProteinComplex(kb, NaN, {'ComplexAA95'}, {'Macromolecular complex AA95'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA96 = ProteinComplex(kb, NaN, {'ComplexAA96'}, {'Macromolecular complex AA96'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA97 = ProteinComplex(kb, NaN, {'ComplexAA97'}, {'Macromolecular complex AA97'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA98 = ProteinComplex(kb, NaN, {'ComplexAA98'}, {'Macromolecular complex AA98'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA99 = ProteinComplex(kb, NaN, {'ComplexAA99'}, {'Macromolecular complex AA99'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			cpxAA100 = ProteinComplex(kb, NaN, {'ComplexAA100'}, {'Macromolecular complex AA100'}, 10, {'dsDNA'}, {'dsDNA'}, {''}, {''}, {''}, {''}, {''}, {''}, {''}, {''});
			genome = kb.genome;
			promoter = 'TATCTACAGAGGCCCGCTTTAAGCATCCACGACCCTACTCACTTCAAAGTGGAACCACCCGTCGACGTGTGTCTTAACCCCTGCCGCTGCAAGGTGTGAG';
			lacI = 'ATGAAACCAGTAACGTTATACGATGTCGCAGAGTATGCCGGTGTCTCTTATCAGACCGTTTCCCGCGTGGTGAACCAGGCCAGCCACGTTTCTGCGAAAACGCGGGAAAAAGTGGAAGCGGCGATGGCGGAGCTGAATTACATTCCCAACCGCGTGGCACAACAACTGGCGGGCAAACAGTCGTTGCTGATTGGCGTTGCCACCTCCAGTCTGGCCCTGCACGCGCCGTCGCAAATTGTCGCGGCGATTAAATCTCGCGCCGATCAACTGGGTGCCAGCGTGGTGGTGTCGATGGTAGAACGAAGCGGCGTCGAAGCCTGTAAAGCGGCGGTGCACAATCTTCTCGCGCAACGCGTCAGTGGGCTGATCATTAACTATCCGCTGGATGACCAGGATGCCATTGCTGTGGAAGCTGCCTGCACTAATGTTCCGGCGTTATTTCTTGATGTCTCTGACCAGACACCCATCAACAGTATTATTTTCTCCCATGAAGACGGTACGCGACTGGGCGTGGAGCATCTGGTCGCATTGGGTCACCAGCAAATCGCGCTGTTAGCGGGCCCATTAAGTTCTGTCTCGGCGCGTCTGCGTCTGGCTGGCTGGCATAAATATCTCACTCGCAATCAAATTCAGCCGATAGCGGAACGGGAAGGCGACTGGAGTGCCATGTCCGGTTTTCAACAAACCATGCAAATGCTGAATGAGGGCATCGTTCCCACTGCGATGCTGGTTGCCAACGATCAGATGGCGCTGGGCGCAATGCGCGCCATTACCGAGTCCGGGCTGCGCGTTGGTGCGGATATCTCGGTAGTGGGATACGACGATACCGAAGACAGCTCATGTTATATCCCGCCGTCAACCACCATCAAACAGGATTTTCGCCTGCTGGGGCAAACCAGCGTGGACCGCTTGCTGCAACTCTCTCAGGGCCAGGCGGTGAAGGGCAATCAGCTGTTGCCCGTCTCACTGGTGAAAAGAAAAACCACCCTGGCGCCCAATACGCAAACCGCCTCTCCCCGCGCGTTGGCCGATTCATTAATGCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCGGGCAGTAA';
			genome.sequence = [genome.sequence repmat([promoter lacI], [1, 100])];
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

			geneA51.genome = genome;
			geneA51.transcriptionUnits = tuA51;
			geneA51.proteinMonomers = monA51;
			geneA51.compartment = cytosol;

			geneA52.genome = genome;
			geneA52.transcriptionUnits = tuA52;
			geneA52.proteinMonomers = monA52;
			geneA52.compartment = cytosol;

			geneA53.genome = genome;
			geneA53.transcriptionUnits = tuA53;
			geneA53.proteinMonomers = monA53;
			geneA53.compartment = cytosol;

			geneA54.genome = genome;
			geneA54.transcriptionUnits = tuA54;
			geneA54.proteinMonomers = monA54;
			geneA54.compartment = cytosol;

			geneA55.genome = genome;
			geneA55.transcriptionUnits = tuA55;
			geneA55.proteinMonomers = monA55;
			geneA55.compartment = cytosol;

			geneA56.genome = genome;
			geneA56.transcriptionUnits = tuA56;
			geneA56.proteinMonomers = monA56;
			geneA56.compartment = cytosol;

			geneA57.genome = genome;
			geneA57.transcriptionUnits = tuA57;
			geneA57.proteinMonomers = monA57;
			geneA57.compartment = cytosol;

			geneA58.genome = genome;
			geneA58.transcriptionUnits = tuA58;
			geneA58.proteinMonomers = monA58;
			geneA58.compartment = cytosol;

			geneA59.genome = genome;
			geneA59.transcriptionUnits = tuA59;
			geneA59.proteinMonomers = monA59;
			geneA59.compartment = cytosol;

			geneA60.genome = genome;
			geneA60.transcriptionUnits = tuA60;
			geneA60.proteinMonomers = monA60;
			geneA60.compartment = cytosol;

			geneA61.genome = genome;
			geneA61.transcriptionUnits = tuA61;
			geneA61.proteinMonomers = monA61;
			geneA61.compartment = cytosol;

			geneA62.genome = genome;
			geneA62.transcriptionUnits = tuA62;
			geneA62.proteinMonomers = monA62;
			geneA62.compartment = cytosol;

			geneA63.genome = genome;
			geneA63.transcriptionUnits = tuA63;
			geneA63.proteinMonomers = monA63;
			geneA63.compartment = cytosol;

			geneA64.genome = genome;
			geneA64.transcriptionUnits = tuA64;
			geneA64.proteinMonomers = monA64;
			geneA64.compartment = cytosol;

			geneA65.genome = genome;
			geneA65.transcriptionUnits = tuA65;
			geneA65.proteinMonomers = monA65;
			geneA65.compartment = cytosol;

			geneA66.genome = genome;
			geneA66.transcriptionUnits = tuA66;
			geneA66.proteinMonomers = monA66;
			geneA66.compartment = cytosol;

			geneA67.genome = genome;
			geneA67.transcriptionUnits = tuA67;
			geneA67.proteinMonomers = monA67;
			geneA67.compartment = cytosol;

			geneA68.genome = genome;
			geneA68.transcriptionUnits = tuA68;
			geneA68.proteinMonomers = monA68;
			geneA68.compartment = cytosol;

			geneA69.genome = genome;
			geneA69.transcriptionUnits = tuA69;
			geneA69.proteinMonomers = monA69;
			geneA69.compartment = cytosol;

			geneA70.genome = genome;
			geneA70.transcriptionUnits = tuA70;
			geneA70.proteinMonomers = monA70;
			geneA70.compartment = cytosol;

			geneA71.genome = genome;
			geneA71.transcriptionUnits = tuA71;
			geneA71.proteinMonomers = monA71;
			geneA71.compartment = cytosol;

			geneA72.genome = genome;
			geneA72.transcriptionUnits = tuA72;
			geneA72.proteinMonomers = monA72;
			geneA72.compartment = cytosol;

			geneA73.genome = genome;
			geneA73.transcriptionUnits = tuA73;
			geneA73.proteinMonomers = monA73;
			geneA73.compartment = cytosol;

			geneA74.genome = genome;
			geneA74.transcriptionUnits = tuA74;
			geneA74.proteinMonomers = monA74;
			geneA74.compartment = cytosol;

			geneA75.genome = genome;
			geneA75.transcriptionUnits = tuA75;
			geneA75.proteinMonomers = monA75;
			geneA75.compartment = cytosol;

			geneA76.genome = genome;
			geneA76.transcriptionUnits = tuA76;
			geneA76.proteinMonomers = monA76;
			geneA76.compartment = cytosol;

			geneA77.genome = genome;
			geneA77.transcriptionUnits = tuA77;
			geneA77.proteinMonomers = monA77;
			geneA77.compartment = cytosol;

			geneA78.genome = genome;
			geneA78.transcriptionUnits = tuA78;
			geneA78.proteinMonomers = monA78;
			geneA78.compartment = cytosol;

			geneA79.genome = genome;
			geneA79.transcriptionUnits = tuA79;
			geneA79.proteinMonomers = monA79;
			geneA79.compartment = cytosol;

			geneA80.genome = genome;
			geneA80.transcriptionUnits = tuA80;
			geneA80.proteinMonomers = monA80;
			geneA80.compartment = cytosol;

			geneA81.genome = genome;
			geneA81.transcriptionUnits = tuA81;
			geneA81.proteinMonomers = monA81;
			geneA81.compartment = cytosol;

			geneA82.genome = genome;
			geneA82.transcriptionUnits = tuA82;
			geneA82.proteinMonomers = monA82;
			geneA82.compartment = cytosol;

			geneA83.genome = genome;
			geneA83.transcriptionUnits = tuA83;
			geneA83.proteinMonomers = monA83;
			geneA83.compartment = cytosol;

			geneA84.genome = genome;
			geneA84.transcriptionUnits = tuA84;
			geneA84.proteinMonomers = monA84;
			geneA84.compartment = cytosol;

			geneA85.genome = genome;
			geneA85.transcriptionUnits = tuA85;
			geneA85.proteinMonomers = monA85;
			geneA85.compartment = cytosol;

			geneA86.genome = genome;
			geneA86.transcriptionUnits = tuA86;
			geneA86.proteinMonomers = monA86;
			geneA86.compartment = cytosol;

			geneA87.genome = genome;
			geneA87.transcriptionUnits = tuA87;
			geneA87.proteinMonomers = monA87;
			geneA87.compartment = cytosol;

			geneA88.genome = genome;
			geneA88.transcriptionUnits = tuA88;
			geneA88.proteinMonomers = monA88;
			geneA88.compartment = cytosol;

			geneA89.genome = genome;
			geneA89.transcriptionUnits = tuA89;
			geneA89.proteinMonomers = monA89;
			geneA89.compartment = cytosol;

			geneA90.genome = genome;
			geneA90.transcriptionUnits = tuA90;
			geneA90.proteinMonomers = monA90;
			geneA90.compartment = cytosol;

			geneA91.genome = genome;
			geneA91.transcriptionUnits = tuA91;
			geneA91.proteinMonomers = monA91;
			geneA91.compartment = cytosol;

			geneA92.genome = genome;
			geneA92.transcriptionUnits = tuA92;
			geneA92.proteinMonomers = monA92;
			geneA92.compartment = cytosol;

			geneA93.genome = genome;
			geneA93.transcriptionUnits = tuA93;
			geneA93.proteinMonomers = monA93;
			geneA93.compartment = cytosol;

			geneA94.genome = genome;
			geneA94.transcriptionUnits = tuA94;
			geneA94.proteinMonomers = monA94;
			geneA94.compartment = cytosol;

			geneA95.genome = genome;
			geneA95.transcriptionUnits = tuA95;
			geneA95.proteinMonomers = monA95;
			geneA95.compartment = cytosol;

			geneA96.genome = genome;
			geneA96.transcriptionUnits = tuA96;
			geneA96.proteinMonomers = monA96;
			geneA96.compartment = cytosol;

			geneA97.genome = genome;
			geneA97.transcriptionUnits = tuA97;
			geneA97.proteinMonomers = monA97;
			geneA97.compartment = cytosol;

			geneA98.genome = genome;
			geneA98.transcriptionUnits = tuA98;
			geneA98.proteinMonomers = monA98;
			geneA98.compartment = cytosol;

			geneA99.genome = genome;
			geneA99.transcriptionUnits = tuA99;
			geneA99.proteinMonomers = monA99;
			geneA99.compartment = cytosol;

			geneA100.genome = genome;
			geneA100.transcriptionUnits = tuA100;
			geneA100.proteinMonomers = monA100;
			geneA100.compartment = cytosol;

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

			tuA51.genome = genome;
			tuA51.genes = geneA51;
			tuA51.geneCompartments = cytosol;
			tuA51.compartment = cytosol;

			tuA52.genome = genome;
			tuA52.genes = geneA52;
			tuA52.geneCompartments = cytosol;
			tuA52.compartment = cytosol;

			tuA53.genome = genome;
			tuA53.genes = geneA53;
			tuA53.geneCompartments = cytosol;
			tuA53.compartment = cytosol;

			tuA54.genome = genome;
			tuA54.genes = geneA54;
			tuA54.geneCompartments = cytosol;
			tuA54.compartment = cytosol;

			tuA55.genome = genome;
			tuA55.genes = geneA55;
			tuA55.geneCompartments = cytosol;
			tuA55.compartment = cytosol;

			tuA56.genome = genome;
			tuA56.genes = geneA56;
			tuA56.geneCompartments = cytosol;
			tuA56.compartment = cytosol;

			tuA57.genome = genome;
			tuA57.genes = geneA57;
			tuA57.geneCompartments = cytosol;
			tuA57.compartment = cytosol;

			tuA58.genome = genome;
			tuA58.genes = geneA58;
			tuA58.geneCompartments = cytosol;
			tuA58.compartment = cytosol;

			tuA59.genome = genome;
			tuA59.genes = geneA59;
			tuA59.geneCompartments = cytosol;
			tuA59.compartment = cytosol;

			tuA60.genome = genome;
			tuA60.genes = geneA60;
			tuA60.geneCompartments = cytosol;
			tuA60.compartment = cytosol;

			tuA61.genome = genome;
			tuA61.genes = geneA61;
			tuA61.geneCompartments = cytosol;
			tuA61.compartment = cytosol;

			tuA62.genome = genome;
			tuA62.genes = geneA62;
			tuA62.geneCompartments = cytosol;
			tuA62.compartment = cytosol;

			tuA63.genome = genome;
			tuA63.genes = geneA63;
			tuA63.geneCompartments = cytosol;
			tuA63.compartment = cytosol;

			tuA64.genome = genome;
			tuA64.genes = geneA64;
			tuA64.geneCompartments = cytosol;
			tuA64.compartment = cytosol;

			tuA65.genome = genome;
			tuA65.genes = geneA65;
			tuA65.geneCompartments = cytosol;
			tuA65.compartment = cytosol;

			tuA66.genome = genome;
			tuA66.genes = geneA66;
			tuA66.geneCompartments = cytosol;
			tuA66.compartment = cytosol;

			tuA67.genome = genome;
			tuA67.genes = geneA67;
			tuA67.geneCompartments = cytosol;
			tuA67.compartment = cytosol;

			tuA68.genome = genome;
			tuA68.genes = geneA68;
			tuA68.geneCompartments = cytosol;
			tuA68.compartment = cytosol;

			tuA69.genome = genome;
			tuA69.genes = geneA69;
			tuA69.geneCompartments = cytosol;
			tuA69.compartment = cytosol;

			tuA70.genome = genome;
			tuA70.genes = geneA70;
			tuA70.geneCompartments = cytosol;
			tuA70.compartment = cytosol;

			tuA71.genome = genome;
			tuA71.genes = geneA71;
			tuA71.geneCompartments = cytosol;
			tuA71.compartment = cytosol;

			tuA72.genome = genome;
			tuA72.genes = geneA72;
			tuA72.geneCompartments = cytosol;
			tuA72.compartment = cytosol;

			tuA73.genome = genome;
			tuA73.genes = geneA73;
			tuA73.geneCompartments = cytosol;
			tuA73.compartment = cytosol;

			tuA74.genome = genome;
			tuA74.genes = geneA74;
			tuA74.geneCompartments = cytosol;
			tuA74.compartment = cytosol;

			tuA75.genome = genome;
			tuA75.genes = geneA75;
			tuA75.geneCompartments = cytosol;
			tuA75.compartment = cytosol;

			tuA76.genome = genome;
			tuA76.genes = geneA76;
			tuA76.geneCompartments = cytosol;
			tuA76.compartment = cytosol;

			tuA77.genome = genome;
			tuA77.genes = geneA77;
			tuA77.geneCompartments = cytosol;
			tuA77.compartment = cytosol;

			tuA78.genome = genome;
			tuA78.genes = geneA78;
			tuA78.geneCompartments = cytosol;
			tuA78.compartment = cytosol;

			tuA79.genome = genome;
			tuA79.genes = geneA79;
			tuA79.geneCompartments = cytosol;
			tuA79.compartment = cytosol;

			tuA80.genome = genome;
			tuA80.genes = geneA80;
			tuA80.geneCompartments = cytosol;
			tuA80.compartment = cytosol;

			tuA81.genome = genome;
			tuA81.genes = geneA81;
			tuA81.geneCompartments = cytosol;
			tuA81.compartment = cytosol;

			tuA82.genome = genome;
			tuA82.genes = geneA82;
			tuA82.geneCompartments = cytosol;
			tuA82.compartment = cytosol;

			tuA83.genome = genome;
			tuA83.genes = geneA83;
			tuA83.geneCompartments = cytosol;
			tuA83.compartment = cytosol;

			tuA84.genome = genome;
			tuA84.genes = geneA84;
			tuA84.geneCompartments = cytosol;
			tuA84.compartment = cytosol;

			tuA85.genome = genome;
			tuA85.genes = geneA85;
			tuA85.geneCompartments = cytosol;
			tuA85.compartment = cytosol;

			tuA86.genome = genome;
			tuA86.genes = geneA86;
			tuA86.geneCompartments = cytosol;
			tuA86.compartment = cytosol;

			tuA87.genome = genome;
			tuA87.genes = geneA87;
			tuA87.geneCompartments = cytosol;
			tuA87.compartment = cytosol;

			tuA88.genome = genome;
			tuA88.genes = geneA88;
			tuA88.geneCompartments = cytosol;
			tuA88.compartment = cytosol;

			tuA89.genome = genome;
			tuA89.genes = geneA89;
			tuA89.geneCompartments = cytosol;
			tuA89.compartment = cytosol;

			tuA90.genome = genome;
			tuA90.genes = geneA90;
			tuA90.geneCompartments = cytosol;
			tuA90.compartment = cytosol;

			tuA91.genome = genome;
			tuA91.genes = geneA91;
			tuA91.geneCompartments = cytosol;
			tuA91.compartment = cytosol;

			tuA92.genome = genome;
			tuA92.genes = geneA92;
			tuA92.geneCompartments = cytosol;
			tuA92.compartment = cytosol;

			tuA93.genome = genome;
			tuA93.genes = geneA93;
			tuA93.geneCompartments = cytosol;
			tuA93.compartment = cytosol;

			tuA94.genome = genome;
			tuA94.genes = geneA94;
			tuA94.geneCompartments = cytosol;
			tuA94.compartment = cytosol;

			tuA95.genome = genome;
			tuA95.genes = geneA95;
			tuA95.geneCompartments = cytosol;
			tuA95.compartment = cytosol;

			tuA96.genome = genome;
			tuA96.genes = geneA96;
			tuA96.geneCompartments = cytosol;
			tuA96.compartment = cytosol;

			tuA97.genome = genome;
			tuA97.genes = geneA97;
			tuA97.geneCompartments = cytosol;
			tuA97.compartment = cytosol;

			tuA98.genome = genome;
			tuA98.genes = geneA98;
			tuA98.geneCompartments = cytosol;
			tuA98.compartment = cytosol;

			tuA99.genome = genome;
			tuA99.genes = geneA99;
			tuA99.geneCompartments = cytosol;
			tuA99.compartment = cytosol;

			tuA100.genome = genome;
			tuA100.genes = geneA100;
			tuA100.geneCompartments = cytosol;
			tuA100.compartment = cytosol;

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

			monA51.gene = geneA51;
			monA51.geneCompartments = cytosol;
			monA51.compartment = cytosol;

			monA52.gene = geneA52;
			monA52.geneCompartments = cytosol;
			monA52.compartment = cytosol;

			monA53.gene = geneA53;
			monA53.geneCompartments = cytosol;
			monA53.compartment = cytosol;

			monA54.gene = geneA54;
			monA54.geneCompartments = cytosol;
			monA54.compartment = cytosol;

			monA55.gene = geneA55;
			monA55.geneCompartments = cytosol;
			monA55.compartment = cytosol;

			monA56.gene = geneA56;
			monA56.geneCompartments = cytosol;
			monA56.compartment = cytosol;

			monA57.gene = geneA57;
			monA57.geneCompartments = cytosol;
			monA57.compartment = cytosol;

			monA58.gene = geneA58;
			monA58.geneCompartments = cytosol;
			monA58.compartment = cytosol;

			monA59.gene = geneA59;
			monA59.geneCompartments = cytosol;
			monA59.compartment = cytosol;

			monA60.gene = geneA60;
			monA60.geneCompartments = cytosol;
			monA60.compartment = cytosol;

			monA61.gene = geneA61;
			monA61.geneCompartments = cytosol;
			monA61.compartment = cytosol;

			monA62.gene = geneA62;
			monA62.geneCompartments = cytosol;
			monA62.compartment = cytosol;

			monA63.gene = geneA63;
			monA63.geneCompartments = cytosol;
			monA63.compartment = cytosol;

			monA64.gene = geneA64;
			monA64.geneCompartments = cytosol;
			monA64.compartment = cytosol;

			monA65.gene = geneA65;
			monA65.geneCompartments = cytosol;
			monA65.compartment = cytosol;

			monA66.gene = geneA66;
			monA66.geneCompartments = cytosol;
			monA66.compartment = cytosol;

			monA67.gene = geneA67;
			monA67.geneCompartments = cytosol;
			monA67.compartment = cytosol;

			monA68.gene = geneA68;
			monA68.geneCompartments = cytosol;
			monA68.compartment = cytosol;

			monA69.gene = geneA69;
			monA69.geneCompartments = cytosol;
			monA69.compartment = cytosol;

			monA70.gene = geneA70;
			monA70.geneCompartments = cytosol;
			monA70.compartment = cytosol;

			monA71.gene = geneA71;
			monA71.geneCompartments = cytosol;
			monA71.compartment = cytosol;

			monA72.gene = geneA72;
			monA72.geneCompartments = cytosol;
			monA72.compartment = cytosol;

			monA73.gene = geneA73;
			monA73.geneCompartments = cytosol;
			monA73.compartment = cytosol;

			monA74.gene = geneA74;
			monA74.geneCompartments = cytosol;
			monA74.compartment = cytosol;

			monA75.gene = geneA75;
			monA75.geneCompartments = cytosol;
			monA75.compartment = cytosol;

			monA76.gene = geneA76;
			monA76.geneCompartments = cytosol;
			monA76.compartment = cytosol;

			monA77.gene = geneA77;
			monA77.geneCompartments = cytosol;
			monA77.compartment = cytosol;

			monA78.gene = geneA78;
			monA78.geneCompartments = cytosol;
			monA78.compartment = cytosol;

			monA79.gene = geneA79;
			monA79.geneCompartments = cytosol;
			monA79.compartment = cytosol;

			monA80.gene = geneA80;
			monA80.geneCompartments = cytosol;
			monA80.compartment = cytosol;

			monA81.gene = geneA81;
			monA81.geneCompartments = cytosol;
			monA81.compartment = cytosol;

			monA82.gene = geneA82;
			monA82.geneCompartments = cytosol;
			monA82.compartment = cytosol;

			monA83.gene = geneA83;
			monA83.geneCompartments = cytosol;
			monA83.compartment = cytosol;

			monA84.gene = geneA84;
			monA84.geneCompartments = cytosol;
			monA84.compartment = cytosol;

			monA85.gene = geneA85;
			monA85.geneCompartments = cytosol;
			monA85.compartment = cytosol;

			monA86.gene = geneA86;
			monA86.geneCompartments = cytosol;
			monA86.compartment = cytosol;

			monA87.gene = geneA87;
			monA87.geneCompartments = cytosol;
			monA87.compartment = cytosol;

			monA88.gene = geneA88;
			monA88.geneCompartments = cytosol;
			monA88.compartment = cytosol;

			monA89.gene = geneA89;
			monA89.geneCompartments = cytosol;
			monA89.compartment = cytosol;

			monA90.gene = geneA90;
			monA90.geneCompartments = cytosol;
			monA90.compartment = cytosol;

			monA91.gene = geneA91;
			monA91.geneCompartments = cytosol;
			monA91.compartment = cytosol;

			monA92.gene = geneA92;
			monA92.geneCompartments = cytosol;
			monA92.compartment = cytosol;

			monA93.gene = geneA93;
			monA93.geneCompartments = cytosol;
			monA93.compartment = cytosol;

			monA94.gene = geneA94;
			monA94.geneCompartments = cytosol;
			monA94.compartment = cytosol;

			monA95.gene = geneA95;
			monA95.geneCompartments = cytosol;
			monA95.compartment = cytosol;

			monA96.gene = geneA96;
			monA96.geneCompartments = cytosol;
			monA96.compartment = cytosol;

			monA97.gene = geneA97;
			monA97.geneCompartments = cytosol;
			monA97.compartment = cytosol;

			monA98.gene = geneA98;
			monA98.geneCompartments = cytosol;
			monA98.compartment = cytosol;

			monA99.gene = geneA99;
			monA99.geneCompartments = cytosol;
			monA99.compartment = cytosol;

			monA100.gene = geneA100;
			monA100.geneCompartments = cytosol;
			monA100.compartment = cytosol;

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

			cpxAA51.proteinMonomers = monA51;
			cpxAA51.proteinMonomerCompartments = cytosol;
			cpxAA51.proteinMonomerCoefficients = 4;
			cpxAA51.compartment = cytosol;
			cpxAA51.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA52.proteinMonomers = monA52;
			cpxAA52.proteinMonomerCompartments = cytosol;
			cpxAA52.proteinMonomerCoefficients = 4;
			cpxAA52.compartment = cytosol;
			cpxAA52.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA53.proteinMonomers = monA53;
			cpxAA53.proteinMonomerCompartments = cytosol;
			cpxAA53.proteinMonomerCoefficients = 4;
			cpxAA53.compartment = cytosol;
			cpxAA53.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA54.proteinMonomers = monA54;
			cpxAA54.proteinMonomerCompartments = cytosol;
			cpxAA54.proteinMonomerCoefficients = 4;
			cpxAA54.compartment = cytosol;
			cpxAA54.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA55.proteinMonomers = monA55;
			cpxAA55.proteinMonomerCompartments = cytosol;
			cpxAA55.proteinMonomerCoefficients = 4;
			cpxAA55.compartment = cytosol;
			cpxAA55.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA56.proteinMonomers = monA56;
			cpxAA56.proteinMonomerCompartments = cytosol;
			cpxAA56.proteinMonomerCoefficients = 4;
			cpxAA56.compartment = cytosol;
			cpxAA56.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA57.proteinMonomers = monA57;
			cpxAA57.proteinMonomerCompartments = cytosol;
			cpxAA57.proteinMonomerCoefficients = 4;
			cpxAA57.compartment = cytosol;
			cpxAA57.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA58.proteinMonomers = monA58;
			cpxAA58.proteinMonomerCompartments = cytosol;
			cpxAA58.proteinMonomerCoefficients = 4;
			cpxAA58.compartment = cytosol;
			cpxAA58.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA59.proteinMonomers = monA59;
			cpxAA59.proteinMonomerCompartments = cytosol;
			cpxAA59.proteinMonomerCoefficients = 4;
			cpxAA59.compartment = cytosol;
			cpxAA59.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA60.proteinMonomers = monA60;
			cpxAA60.proteinMonomerCompartments = cytosol;
			cpxAA60.proteinMonomerCoefficients = 4;
			cpxAA60.compartment = cytosol;
			cpxAA60.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA61.proteinMonomers = monA61;
			cpxAA61.proteinMonomerCompartments = cytosol;
			cpxAA61.proteinMonomerCoefficients = 4;
			cpxAA61.compartment = cytosol;
			cpxAA61.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA62.proteinMonomers = monA62;
			cpxAA62.proteinMonomerCompartments = cytosol;
			cpxAA62.proteinMonomerCoefficients = 4;
			cpxAA62.compartment = cytosol;
			cpxAA62.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA63.proteinMonomers = monA63;
			cpxAA63.proteinMonomerCompartments = cytosol;
			cpxAA63.proteinMonomerCoefficients = 4;
			cpxAA63.compartment = cytosol;
			cpxAA63.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA64.proteinMonomers = monA64;
			cpxAA64.proteinMonomerCompartments = cytosol;
			cpxAA64.proteinMonomerCoefficients = 4;
			cpxAA64.compartment = cytosol;
			cpxAA64.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA65.proteinMonomers = monA65;
			cpxAA65.proteinMonomerCompartments = cytosol;
			cpxAA65.proteinMonomerCoefficients = 4;
			cpxAA65.compartment = cytosol;
			cpxAA65.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA66.proteinMonomers = monA66;
			cpxAA66.proteinMonomerCompartments = cytosol;
			cpxAA66.proteinMonomerCoefficients = 4;
			cpxAA66.compartment = cytosol;
			cpxAA66.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA67.proteinMonomers = monA67;
			cpxAA67.proteinMonomerCompartments = cytosol;
			cpxAA67.proteinMonomerCoefficients = 4;
			cpxAA67.compartment = cytosol;
			cpxAA67.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA68.proteinMonomers = monA68;
			cpxAA68.proteinMonomerCompartments = cytosol;
			cpxAA68.proteinMonomerCoefficients = 4;
			cpxAA68.compartment = cytosol;
			cpxAA68.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA69.proteinMonomers = monA69;
			cpxAA69.proteinMonomerCompartments = cytosol;
			cpxAA69.proteinMonomerCoefficients = 4;
			cpxAA69.compartment = cytosol;
			cpxAA69.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA70.proteinMonomers = monA70;
			cpxAA70.proteinMonomerCompartments = cytosol;
			cpxAA70.proteinMonomerCoefficients = 4;
			cpxAA70.compartment = cytosol;
			cpxAA70.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA71.proteinMonomers = monA71;
			cpxAA71.proteinMonomerCompartments = cytosol;
			cpxAA71.proteinMonomerCoefficients = 4;
			cpxAA71.compartment = cytosol;
			cpxAA71.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA72.proteinMonomers = monA72;
			cpxAA72.proteinMonomerCompartments = cytosol;
			cpxAA72.proteinMonomerCoefficients = 4;
			cpxAA72.compartment = cytosol;
			cpxAA72.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA73.proteinMonomers = monA73;
			cpxAA73.proteinMonomerCompartments = cytosol;
			cpxAA73.proteinMonomerCoefficients = 4;
			cpxAA73.compartment = cytosol;
			cpxAA73.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA74.proteinMonomers = monA74;
			cpxAA74.proteinMonomerCompartments = cytosol;
			cpxAA74.proteinMonomerCoefficients = 4;
			cpxAA74.compartment = cytosol;
			cpxAA74.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA75.proteinMonomers = monA75;
			cpxAA75.proteinMonomerCompartments = cytosol;
			cpxAA75.proteinMonomerCoefficients = 4;
			cpxAA75.compartment = cytosol;
			cpxAA75.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA76.proteinMonomers = monA76;
			cpxAA76.proteinMonomerCompartments = cytosol;
			cpxAA76.proteinMonomerCoefficients = 4;
			cpxAA76.compartment = cytosol;
			cpxAA76.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA77.proteinMonomers = monA77;
			cpxAA77.proteinMonomerCompartments = cytosol;
			cpxAA77.proteinMonomerCoefficients = 4;
			cpxAA77.compartment = cytosol;
			cpxAA77.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA78.proteinMonomers = monA78;
			cpxAA78.proteinMonomerCompartments = cytosol;
			cpxAA78.proteinMonomerCoefficients = 4;
			cpxAA78.compartment = cytosol;
			cpxAA78.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA79.proteinMonomers = monA79;
			cpxAA79.proteinMonomerCompartments = cytosol;
			cpxAA79.proteinMonomerCoefficients = 4;
			cpxAA79.compartment = cytosol;
			cpxAA79.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA80.proteinMonomers = monA80;
			cpxAA80.proteinMonomerCompartments = cytosol;
			cpxAA80.proteinMonomerCoefficients = 4;
			cpxAA80.compartment = cytosol;
			cpxAA80.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA81.proteinMonomers = monA81;
			cpxAA81.proteinMonomerCompartments = cytosol;
			cpxAA81.proteinMonomerCoefficients = 4;
			cpxAA81.compartment = cytosol;
			cpxAA81.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA82.proteinMonomers = monA82;
			cpxAA82.proteinMonomerCompartments = cytosol;
			cpxAA82.proteinMonomerCoefficients = 4;
			cpxAA82.compartment = cytosol;
			cpxAA82.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA83.proteinMonomers = monA83;
			cpxAA83.proteinMonomerCompartments = cytosol;
			cpxAA83.proteinMonomerCoefficients = 4;
			cpxAA83.compartment = cytosol;
			cpxAA83.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA84.proteinMonomers = monA84;
			cpxAA84.proteinMonomerCompartments = cytosol;
			cpxAA84.proteinMonomerCoefficients = 4;
			cpxAA84.compartment = cytosol;
			cpxAA84.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA85.proteinMonomers = monA85;
			cpxAA85.proteinMonomerCompartments = cytosol;
			cpxAA85.proteinMonomerCoefficients = 4;
			cpxAA85.compartment = cytosol;
			cpxAA85.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA86.proteinMonomers = monA86;
			cpxAA86.proteinMonomerCompartments = cytosol;
			cpxAA86.proteinMonomerCoefficients = 4;
			cpxAA86.compartment = cytosol;
			cpxAA86.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA87.proteinMonomers = monA87;
			cpxAA87.proteinMonomerCompartments = cytosol;
			cpxAA87.proteinMonomerCoefficients = 4;
			cpxAA87.compartment = cytosol;
			cpxAA87.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA88.proteinMonomers = monA88;
			cpxAA88.proteinMonomerCompartments = cytosol;
			cpxAA88.proteinMonomerCoefficients = 4;
			cpxAA88.compartment = cytosol;
			cpxAA88.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA89.proteinMonomers = monA89;
			cpxAA89.proteinMonomerCompartments = cytosol;
			cpxAA89.proteinMonomerCoefficients = 4;
			cpxAA89.compartment = cytosol;
			cpxAA89.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA90.proteinMonomers = monA90;
			cpxAA90.proteinMonomerCompartments = cytosol;
			cpxAA90.proteinMonomerCoefficients = 4;
			cpxAA90.compartment = cytosol;
			cpxAA90.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA91.proteinMonomers = monA91;
			cpxAA91.proteinMonomerCompartments = cytosol;
			cpxAA91.proteinMonomerCoefficients = 4;
			cpxAA91.compartment = cytosol;
			cpxAA91.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA92.proteinMonomers = monA92;
			cpxAA92.proteinMonomerCompartments = cytosol;
			cpxAA92.proteinMonomerCoefficients = 4;
			cpxAA92.compartment = cytosol;
			cpxAA92.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA93.proteinMonomers = monA93;
			cpxAA93.proteinMonomerCompartments = cytosol;
			cpxAA93.proteinMonomerCoefficients = 4;
			cpxAA93.compartment = cytosol;
			cpxAA93.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA94.proteinMonomers = monA94;
			cpxAA94.proteinMonomerCompartments = cytosol;
			cpxAA94.proteinMonomerCoefficients = 4;
			cpxAA94.compartment = cytosol;
			cpxAA94.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA95.proteinMonomers = monA95;
			cpxAA95.proteinMonomerCompartments = cytosol;
			cpxAA95.proteinMonomerCoefficients = 4;
			cpxAA95.compartment = cytosol;
			cpxAA95.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA96.proteinMonomers = monA96;
			cpxAA96.proteinMonomerCompartments = cytosol;
			cpxAA96.proteinMonomerCoefficients = 4;
			cpxAA96.compartment = cytosol;
			cpxAA96.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA97.proteinMonomers = monA97;
			cpxAA97.proteinMonomerCompartments = cytosol;
			cpxAA97.proteinMonomerCoefficients = 4;
			cpxAA97.compartment = cytosol;
			cpxAA97.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA98.proteinMonomers = monA98;
			cpxAA98.proteinMonomerCompartments = cytosol;
			cpxAA98.proteinMonomerCoefficients = 4;
			cpxAA98.compartment = cytosol;
			cpxAA98.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA99.proteinMonomers = monA99;
			cpxAA99.proteinMonomerCompartments = cytosol;
			cpxAA99.proteinMonomerCoefficients = 4;
			cpxAA99.compartment = cytosol;
			cpxAA99.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			cpxAA100.proteinMonomers = monA100;
			cpxAA100.proteinMonomerCompartments = cytosol;
			cpxAA100.proteinMonomerCoefficients = 4;
			cpxAA100.compartment = cytosol;
			cpxAA100.complexFormationProcess = kb.processes.findobj('wholeCellModelID', 'Process_MacromolecularComplexation');

			genome.genes = [genome.genes; geneA1; geneA2; geneA3; geneA4; geneA5; geneA6; geneA7; geneA8; geneA9; geneA10; geneA11; geneA12; geneA13; geneA14; geneA15; geneA16; geneA17; geneA18; geneA19; geneA20; geneA21; geneA22; geneA23; geneA24; geneA25; geneA26; geneA27; geneA28; geneA29; geneA30; geneA31; geneA32; geneA33; geneA34; geneA35; geneA36; geneA37; geneA38; geneA39; geneA40; geneA41; geneA42; geneA43; geneA44; geneA45; geneA46; geneA47; geneA48; geneA49; geneA50; geneA51; geneA52; geneA53; geneA54; geneA55; geneA56; geneA57; geneA58; geneA59; geneA60; geneA61; geneA62; geneA63; geneA64; geneA65; geneA66; geneA67; geneA68; geneA69; geneA70; geneA71; geneA72; geneA73; geneA74; geneA75; geneA76; geneA77; geneA78; geneA79; geneA80; geneA81; geneA82; geneA83; geneA84; geneA85; geneA86; geneA87; geneA88; geneA89; geneA90; geneA91; geneA92; geneA93; geneA94; geneA95; geneA96; geneA97; geneA98; geneA99; geneA100];
			genome.transcriptionUnits = [genome.transcriptionUnits; tuA1; tuA2; tuA3; tuA4; tuA5; tuA6; tuA7; tuA8; tuA9; tuA10; tuA11; tuA12; tuA13; tuA14; tuA15; tuA16; tuA17; tuA18; tuA19; tuA20; tuA21; tuA22; tuA23; tuA24; tuA25; tuA26; tuA27; tuA28; tuA29; tuA30; tuA31; tuA32; tuA33; tuA34; tuA35; tuA36; tuA37; tuA38; tuA39; tuA40; tuA41; tuA42; tuA43; tuA44; tuA45; tuA46; tuA47; tuA48; tuA49; tuA50; tuA51; tuA52; tuA53; tuA54; tuA55; tuA56; tuA57; tuA58; tuA59; tuA60; tuA61; tuA62; tuA63; tuA64; tuA65; tuA66; tuA67; tuA68; tuA69; tuA70; tuA71; tuA72; tuA73; tuA74; tuA75; tuA76; tuA77; tuA78; tuA79; tuA80; tuA81; tuA82; tuA83; tuA84; tuA85; tuA86; tuA87; tuA88; tuA89; tuA90; tuA91; tuA92; tuA93; tuA94; tuA95; tuA96; tuA97; tuA98; tuA99; tuA100];
			kb.genes = genome.genes;
			kb.transcriptionUnits = genome.transcriptionUnits;
			kb.proteinMonomers = [kb.proteinMonomers; monA1; monA2; monA3; monA4; monA5; monA6; monA7; monA8; monA9; monA10; monA11; monA12; monA13; monA14; monA15; monA16; monA17; monA18; monA19; monA20; monA21; monA22; monA23; monA24; monA25; monA26; monA27; monA28; monA29; monA30; monA31; monA32; monA33; monA34; monA35; monA36; monA37; monA38; monA39; monA40; monA41; monA42; monA43; monA44; monA45; monA46; monA47; monA48; monA49; monA50; monA51; monA52; monA53; monA54; monA55; monA56; monA57; monA58; monA59; monA60; monA61; monA62; monA63; monA64; monA65; monA66; monA67; monA68; monA69; monA70; monA71; monA72; monA73; monA74; monA75; monA76; monA77; monA78; monA79; monA80; monA81; monA82; monA83; monA84; monA85; monA86; monA87; monA88; monA89; monA90; monA91; monA92; monA93; monA94; monA95; monA96; monA97; monA98; monA99; monA100];
			kb.proteinComplexs = [kb.proteinComplexs; cpxAA1; cpxAA2; cpxAA3; cpxAA4; cpxAA5; cpxAA6; cpxAA7; cpxAA8; cpxAA9; cpxAA10; cpxAA11; cpxAA12; cpxAA13; cpxAA14; cpxAA15; cpxAA16; cpxAA17; cpxAA18; cpxAA19; cpxAA20; cpxAA21; cpxAA22; cpxAA23; cpxAA24; cpxAA25; cpxAA26; cpxAA27; cpxAA28; cpxAA29; cpxAA30; cpxAA31; cpxAA32; cpxAA33; cpxAA34; cpxAA35; cpxAA36; cpxAA37; cpxAA38; cpxAA39; cpxAA40; cpxAA41; cpxAA42; cpxAA43; cpxAA44; cpxAA45; cpxAA46; cpxAA47; cpxAA48; cpxAA49; cpxAA50; cpxAA51; cpxAA52; cpxAA53; cpxAA54; cpxAA55; cpxAA56; cpxAA57; cpxAA58; cpxAA59; cpxAA60; cpxAA61; cpxAA62; cpxAA63; cpxAA64; cpxAA65; cpxAA66; cpxAA67; cpxAA68; cpxAA69; cpxAA70; cpxAA71; cpxAA72; cpxAA73; cpxAA74; cpxAA75; cpxAA76; cpxAA77; cpxAA78; cpxAA79; cpxAA80; cpxAA81; cpxAA82; cpxAA83; cpxAA84; cpxAA85; cpxAA86; cpxAA87; cpxAA88; cpxAA89; cpxAA90; cpxAA91; cpxAA92; cpxAA93; cpxAA94; cpxAA95; cpxAA96; cpxAA97; cpxAA98; cpxAA99; cpxAA100];
			this.modifyNetworkStructure@edu.stanford.covert.cell.sim.runners.SimulationRunner(kb);

		end		function modifyNetworkParameters(~, sim)
			g = sim.gene;
			time = sim.state('Time');
			rna = sim.state('Rna');
			trn = sim.process('Transcription');
			nascentMRNAIndexs = find(any(rna.nascentRNAGeneComposition(g.mRNAIndexs, :), 1));
			[~, modTuIndexs] = ismember({'TuA1', 'TuA2', 'TuA3', 'TuA4', 'TuA5', 'TuA6', 'TuA7', 'TuA8', 'TuA9', 'TuA10', 'TuA11', 'TuA12', 'TuA13', 'TuA14', 'TuA15', 'TuA16', 'TuA17', 'TuA18', 'TuA19', 'TuA20', 'TuA21', 'TuA22', 'TuA23', 'TuA24', 'TuA25', 'TuA26', 'TuA27', 'TuA28', 'TuA29', 'TuA30', 'TuA31', 'TuA32', 'TuA33', 'TuA34', 'TuA35', 'TuA36', 'TuA37', 'TuA38', 'TuA39', 'TuA40', 'TuA41', 'TuA42', 'TuA43', 'TuA44', 'TuA45', 'TuA46', 'TuA47', 'TuA48', 'TuA49', 'TuA50', 'TuA51', 'TuA52', 'TuA53', 'TuA54', 'TuA55', 'TuA56', 'TuA57', 'TuA58', 'TuA59', 'TuA60', 'TuA61', 'TuA62', 'TuA63', 'TuA64', 'TuA65', 'TuA66', 'TuA67', 'TuA68', 'TuA69', 'TuA70', 'TuA71', 'TuA72', 'TuA73', 'TuA74', 'TuA75', 'TuA76', 'TuA77', 'TuA78', 'TuA79', 'TuA80', 'TuA81', 'TuA82', 'TuA83', 'TuA84', 'TuA85', 'TuA86', 'TuA87', 'TuA88', 'TuA89', 'TuA90', 'TuA91', 'TuA92', 'TuA93', 'TuA94', 'TuA95', 'TuA96', 'TuA97', 'TuA98', 'TuA99', 'TuA100'}, rna.wholeCellModelIDs(rna.nascentIndexs));
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