% EmpiricalFormula test cases
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 9/5/2010
classdef EmpiricalFormula_Test < TestCase
    %constructor
    methods
        function this = EmpiricalFormula_Test(name)
            this = this@TestCase(name);
        end
    end

    %tests
    methods
        function testConstructor(~)
            import edu.stanford.covert.util.ConstantUtil;
            import edu.stanford.covert.util.EmpiricalFormula;

            %% Null
            %ex 1
            ef = EmpiricalFormula('');
            assertEqual([0 0], size(ef));

            %% char
            %ex 1
            ef = EmpiricalFormula('H2O');
            assertEqual([1 1], size(ef));
            assertEqual(2, ef.H);
            assertEqual(1, ef.O);
            assertEqual(struct, ef.rGroups);
            assertEqual(2*ConstantUtil.elements.H + 1*ConstantUtil.elements.O, ef.mass);
            assertEqual(2*ConstantUtil.elements.H + 1*ConstantUtil.elements.O, ef.weight);
            assertEqual(2*ConstantUtil.elementNumbers.H + 1*ConstantUtil.elementNumbers.O, ef.number);
            assertEqual('H2O', char(ef));

            %ex 2
            ef = EmpiricalFormula('HCl');
            assertEqual([1 1], size(ef));
            assertEqual(1, ef.H);
            assertEqual(1, ef.Cl);
            assertEqual(struct, ef.rGroups);
            assertEqual(1*ConstantUtil.elements.H + 1*ConstantUtil.elements.Cl, ef.mass);
            assertEqual(1*ConstantUtil.elements.H + 1*ConstantUtil.elements.Cl, ef.weight);
            assertEqual(1*ConstantUtil.elementNumbers.H + 1*ConstantUtil.elementNumbers.Cl, ef.number);
            assertEqual('HCl', char(ef));

            %ex 3
            ef = EmpiricalFormula('H1Cl1');
            assertEqual(EmpiricalFormula('HCl'), ef);

            %ex 4
            ef = EmpiricalFormula('C5H12O2');
            assertEqual([1 1], size(ef));
            assertEqual(5, ef.C);
            assertEqual(12, ef.H);
            assertEqual(2, ef.O);
            assertEqual(struct, ef.rGroups);
            assertEqual(12*ConstantUtil.elements.H + 5*ConstantUtil.elements.C + 2*ConstantUtil.elements.O, ef.mass);
            assertEqual(12*ConstantUtil.elements.H + 5*ConstantUtil.elements.C + 2*ConstantUtil.elements.O, ef.weight);
            assertEqual(12*ConstantUtil.elementNumbers.H + 5*ConstantUtil.elementNumbers.C + 2*ConstantUtil.elementNumbers.O, ef.number);
            assertEqual('H12C5O2', char(ef));

            %ex 5
            ef = EmpiricalFormula('ROH');
            assertEqual([1 1], size(ef));
            assertEqual(1, ef.H);
            assertEqual(1, ef.O);
            assertEqual(struct('R',1), ef.rGroups);
            assertEqual(NaN, ef.mass);
            assertEqual(NaN, ef.weight);
            assertEqual(NaN, ef.number);
            assertEqual('HOR', char(ef));

            %% cell array
            %ex1
            ef = EmpiricalFormula({'H2O2'});
            assertEqual(EmpiricalFormula('H2O2'), ef);

            %ex2
            ef = EmpiricalFormula({'H2O2','H','ROH2'});
            assertEqual([1 3], size(ef));
            assertEqual(EmpiricalFormula('H2O2'), ef(1));
            assertEqual(EmpiricalFormula('H'), ef(2));
            assertEqual(EmpiricalFormula('ROH2'), ef(3));

            %ex 3
            ef = EmpiricalFormula({'H2O2';'H';'ROH2'});
            assertEqual([3 1], size(ef));
            assertEqual(EmpiricalFormula('H2O2'), ef(1));
            assertEqual(EmpiricalFormula('H'), ef(2));
            assertEqual(EmpiricalFormula('ROH2'), ef(3));

            %ex 4
            ef = EmpiricalFormula({'H2O2','H';'ROH2','ROH'});
            assertEqual([2 2], size(ef));
            assertEqual(EmpiricalFormula('H2O2'), ef(1,1));
            assertEqual(EmpiricalFormula('H'), ef(1,2));
            assertEqual(EmpiricalFormula('ROH2'), ef(2,1));
            assertEqual(EmpiricalFormula('ROH'), ef(2,2));

            %% struct
            %ex 1
            ef = EmpiricalFormula(struct('H',2,'O',1));
            assertEqual(EmpiricalFormula('H2O'), ef);

            %ex 2
            ef = EmpiricalFormula(repmat(struct('H',2,'O',1),3,1));
            assertEqual(EmpiricalFormula('H2O'), ef(1));
            assertEqual(repmat(EmpiricalFormula('H2O'),3,1), ef);

            %ex 3
            ef = EmpiricalFormula(repmat(struct('H',2,'O',1),1,3));
            assertEqual(repmat(EmpiricalFormula('H2O'),1,3), ef);

            %ex 4
            ef = EmpiricalFormula(repmat(struct('H',2,'O',1),2,2));
            assertEqual(repmat(EmpiricalFormula('H2O'),2,2), ef);
        end

        function testEquality(~)
            import edu.stanford.covert.util.EmpiricalFormula;

            %isequal
            assertTrue(isequal(EmpiricalFormula(''), EmpiricalFormula('')));
            assertTrue(isequal(EmpiricalFormula('H2O'), EmpiricalFormula('H2O1')));
            assertTrue(isequal(EmpiricalFormula('H2O'), EmpiricalFormula({'H2O'})));
            assertTrue(isequal(EmpiricalFormula('H2O'), EmpiricalFormula(struct('H',2,'O',1))));

            assertFalse(isequal(EmpiricalFormula('H2O2'), 'H2O2'));
            assertFalse(isequal(EmpiricalFormula('H2O2'), {'H2O2'}));
            assertFalse(isequal(EmpiricalFormula('H2O2'), struct('H',2,'O',1)));
            assertFalse(isequal(EmpiricalFormula('H2O2'), EmpiricalFormula('H2O')));
            assertFalse(isequal(EmpiricalFormula('H2O2'), EmpiricalFormula('H2O2R')));
            assertFalse(isequal(EmpiricalFormula('H2O2'), EmpiricalFormula({'H2O2','H2O2'})));

            %eq
            assertTrue('H2O' == EmpiricalFormula('H2O'));
            assertTrue(EmpiricalFormula('H2O') == 'H2O');
            assertTrue({'H2O'} == EmpiricalFormula('H2O'));
            assertTrue(EmpiricalFormula('H2O') == {'H2O'});
            assertTrue(struct('H',2,'O',1) == EmpiricalFormula('H2O'));
            assertTrue(EmpiricalFormula('H2O') == struct('H',2,'O',1));

            assertEqual(true(0,0), EmpiricalFormula('') == EmpiricalFormula(''));
            assertEqual(true, EmpiricalFormula('H2O') == EmpiricalFormula('H2O'));
            assertEqual([true false], EmpiricalFormula('H2O') == EmpiricalFormula({'H2O','H2O2'}));
            assertEqual([false true], EmpiricalFormula({'H2O2','H2O'}) == EmpiricalFormula('H2O'));
            assertEqual([false false], EmpiricalFormula('H2OR') == EmpiricalFormula({'H2O','H2O2'}));
            assertEqual([false;false;true], EmpiricalFormula('H2OR') == EmpiricalFormula({'H2O';'H2O2';'H2OR'}));

            %ne
            assertEqual(false(0,0), EmpiricalFormula('') ~= EmpiricalFormula(''));
            assertEqual(false, EmpiricalFormula('H2O') ~= EmpiricalFormula('H2O'));
            assertEqual([false true], EmpiricalFormula('H2O') ~= EmpiricalFormula({'H2O','H2O2'}));
            assertEqual([true false], EmpiricalFormula({'H2O2','H2O'}) ~= EmpiricalFormula('H2O'));
            assertEqual([true true], EmpiricalFormula('H2OR') ~= EmpiricalFormula({'H2O','H2O2'}));
            assertEqual([true;true;false], EmpiricalFormula('H2OR') ~= EmpiricalFormula({'H2O';'H2O2';'H2OR'}));
        end

        function testAdditionSubtraction(~)
            import edu.stanford.covert.util.EmpiricalFormula;

            %uplus
            assertEqual(EmpiricalFormula('H2O'),+EmpiricalFormula('H2O'));

            %uminus
            assertEqual(EmpiricalFormula(struct('H',-2,'O',-1)),-EmpiricalFormula('H2O'));
            assertEqual(EmpiricalFormula(repmat(struct('H',-2,'O',-2,'R',-3),3,3)),-EmpiricalFormula(repmat(struct('H',2,'O',2,'R',3),3,3)));

            %addition
            assertEqual(EmpiricalFormula(struct('H',4,'O',6,'R',2)),EmpiricalFormula(struct('H',2,'O',4,'R',1))+EmpiricalFormula(struct('H',2,'O',2,'R',1)));
            assertEqual(EmpiricalFormula(struct('H',4,'O',6,'R',2)),EmpiricalFormula(struct('H',2,'O',4))+EmpiricalFormula(struct('H',2,'O',2,'R',2)));
            assertEqual(EmpiricalFormula(struct('H',4,'O',6,'R',2)),EmpiricalFormula(struct('H',2,'O',6,'R',1))+EmpiricalFormula(struct('H',2,'R',1)));
            assertEqual(EmpiricalFormula(struct('H',4,'O',6,'R',2)),EmpiricalFormula(struct('H',3,'O',4,'R',1))+EmpiricalFormula(struct('H',1,'O',2,'R',1)));
            assertEqual(repmat(EmpiricalFormula(struct('H',4,'O',6,'R',2)),3,2),repmat(EmpiricalFormula(struct('H',2,'O',4,'R',1)),3,2)+EmpiricalFormula(struct('H',2,'O',2,'R',1)));

            %subtraction
            assertEqual(EmpiricalFormula(struct('H',0,'O',2,'R',0)), EmpiricalFormula(struct('H',2,'O',4,'R',1))-EmpiricalFormula(struct('H',2,'O',2,'R',1)));
            assertEqual(EmpiricalFormula(struct('H',0,'O',2,'R',-2)),EmpiricalFormula(struct('H',2,'O',4))-EmpiricalFormula(struct('H',2,'O',2,'R',2)));
            assertEqual(EmpiricalFormula(struct('H',0,'O',6,'R',0)),EmpiricalFormula(struct('H',2,'O',6,'R',1))-EmpiricalFormula(struct('H',2,'R',1)));
            assertEqual(EmpiricalFormula(struct('H',2,'O',2,'R',0)),EmpiricalFormula(struct('H',3,'O',4,'R',1))-EmpiricalFormula(struct('H',1,'O',2,'R',1)));
            assertEqual(repmat(EmpiricalFormula(struct('H',0,'O',2,'R',0)),3,2),repmat(EmpiricalFormula(struct('H',2,'O',4,'R',1)),3,2)-EmpiricalFormula(struct('H',2,'O',2,'R',1)));
        end

        function testMultiplicationDivision(~)
            import edu.stanford.covert.util.EmpiricalFormula;

            %multiplication
            assertEqual(repmat(EmpiricalFormula(struct('H',4,'O',6,'R',2)),3,2),2*repmat(EmpiricalFormula(struct('H',2,'O',3,'R',1)),3,2));
            assertEqual(repmat(EmpiricalFormula(struct('H',4,'O',6,'R',2)),3,2),2.*repmat(EmpiricalFormula(struct('H',2,'O',3,'R',1)),3,2));
            assertEqual(repmat(EmpiricalFormula(struct('H',4,'O',6,'R',2)),3,2),repmat(EmpiricalFormula(struct('H',2,'O',3,'R',1)),3,2)*2);
            assertEqual(repmat(EmpiricalFormula(struct('H',4,'O',6,'R',2)),3,2),repmat(EmpiricalFormula(struct('H',2,'O',3,'R',1)),3,2).*2);

            %division
            assertEqual(repmat(EmpiricalFormula(struct('H',4,'O',6,'R',2)),3,2),repmat(EmpiricalFormula(struct('H',2,'O',3,'R',1)),3,2)/0.5);
            assertEqual(repmat(EmpiricalFormula(struct('H',4,'O',6,'R',2)),3,2),repmat(EmpiricalFormula(struct('H',2,'O',3,'R',1)),3,2)./0.5);
            assertEqual(repmat(EmpiricalFormula(struct('H',4,'O',6,'R',2)),3,2),0.5\repmat(EmpiricalFormula(struct('H',2,'O',3,'R',1)),3,2));
            assertEqual(repmat(EmpiricalFormula(struct('H',4,'O',6,'R',2)),3,2),0.5.\repmat(EmpiricalFormula(struct('H',2,'O',3,'R',1)),3,2));
        end

        function testSum(~)
            import edu.stanford.covert.util.EmpiricalFormula;

            assertEqual(EmpiricalFormula(struct('H',12,'O',24,'R',18)),sum(repmat(EmpiricalFormula('H2O4R3'),3,2)));
        end

        function testWeightMassNumber(~)
            import edu.stanford.covert.util.ConstantUtil;
            import edu.stanford.covert.util.EmpiricalFormula;

            %ex 1
            ef = EmpiricalFormula('H2O');
            assertEqual(2*ConstantUtil.elements.H + 1*ConstantUtil.elements.O, ef.mass);
            assertEqual(2*ConstantUtil.elements.H + 1*ConstantUtil.elements.O, ef.weight);
            assertEqual(2*ConstantUtil.elementNumbers.H + 1*ConstantUtil.elementNumbers.O, ef.number);
            assertEqual('H2O', char(ef));

            %ex 2
            ef = EmpiricalFormula('HCl');
            assertEqual(1*ConstantUtil.elements.H + 1*ConstantUtil.elements.Cl, ef.mass);
            assertEqual(1*ConstantUtil.elements.H + 1*ConstantUtil.elements.Cl, ef.weight);
            assertEqual(1*ConstantUtil.elementNumbers.H + 1*ConstantUtil.elementNumbers.Cl, ef.number);
            assertEqual('HCl', char(ef));

            %ex 3
            ef = EmpiricalFormula('C5H12O2');
            assertEqual(12*ConstantUtil.elements.H + 5*ConstantUtil.elements.C + 2*ConstantUtil.elements.O, ef.mass);
            assertEqual(12*ConstantUtil.elements.H + 5*ConstantUtil.elements.C + 2*ConstantUtil.elements.O, ef.weight);
            assertEqual(12*ConstantUtil.elementNumbers.H + 5*ConstantUtil.elementNumbers.C + 2*ConstantUtil.elementNumbers.O, ef.number);
            assertEqual('H12C5O2', char(ef));

            %ex 4
            ef = EmpiricalFormula('ROH');
            assertEqual(NaN, ef.mass);
            assertEqual(NaN, ef.weight);
            assertEqual(NaN, ef.number);
            assertEqual('HOR', char(ef));
        end

        function testTypeCasting(~)
            import edu.stanford.covert.util.EmpiricalFormula;

            assertEqual('H2O4R3',cast(EmpiricalFormula('H2O4R3'),'char'));
            assertEqual(repmat(struct('H',2,'O',4,'R',3),3,2),cast(repmat(EmpiricalFormula('H2O4R3'),3,2),'struct'));
            assertEqual({'H2O4R3' 'H2O4R3';'H2O4R3' 'H2O4R3'; 'H2O4R3' 'H2O4R3'},cast(repmat(EmpiricalFormula('H2O4R3'),3,2),'cell'));

            assertEqual('H2O4R3',char(EmpiricalFormula('H2O4R3')));
            assertEqual(repmat(struct('H',2,'O',4,'R',3),3,2),struct(repmat(EmpiricalFormula('H2O4R3'),3,2)));
            assertEqual({'H2O4R3' 'H2O4R3';'H2O4R3' 'H2O4R3'; 'H2O4R3' 'H2O4R3'},cell(repmat(EmpiricalFormula('H2O4R3'),3,2)));
        end

        function testDisplay(~)
            import edu.stanford.covert.util.EmpiricalFormula;
            disp(EmpiricalFormula('H2O4R3'));
            disp(repmat(EmpiricalFormula('H2O4R3'),3,2));
        end
    end
end