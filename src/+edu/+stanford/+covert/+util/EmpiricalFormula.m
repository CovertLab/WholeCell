% EmpiricalFormula
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/18/2010
classdef EmpiricalFormula
    properties
        H  = 0;
        He = 0;
        Li = 0;
        Be = 0;
        B  = 0;
        C  = 0;
        N  = 0;
        O  = 0;
        F  = 0;
        Ne = 0;
        Na = 0;
        Mg = 0;
        Al = 0;
        Si = 0;
        P  = 0;
        S  = 0;
        Cl = 0;
        Ar = 0;
        K  = 0;
        Ca = 0;
        Sc = 0;
        Ti = 0;
        V  = 0;
        Cr = 0;
        Mn = 0;
        Fe = 0;
        Co = 0;
        Ni = 0;
        Cu = 0;
        Zn = 0;
        Ga = 0;
        Ge = 0;
        As = 0;
        Se = 0;
        Br = 0;
        Kr = 0;
        Rb = 0;
        Sr = 0;
        Y  = 0;
        Zr = 0;
        Nb = 0;
        Mo = 0;
        Tc = 0;
        Ru = 0;
        Rh = 0;
        Pd = 0;
        Ag = 0;
        Cd = 0;
        In = 0;
        Sn = 0;
        Sb = 0;
        Te = 0;
        I  = 0;
        Xe = 0;
        Cs = 0;
        Ba = 0;
        La = 0;
        Ce = 0;
        Pr = 0;
        Nd = 0;
        Pm = 0;
        Sm = 0;
        Eu = 0;
        Gd = 0;
        Tb = 0;
        Dy = 0;
        Ho = 0;
        Er = 0;
        Tm = 0;
        Yb = 0;
        Lu = 0;
        Hf = 0;
        Ta = 0;
        W  = 0;
        Re = 0;
        Os = 0;
        Ir = 0;
        Pt = 0;
        Au = 0;
        Hg = 0;
        Tl = 0;
        Pb = 0;
        Bi = 0;
        Po = 0;
        At = 0;
        Rn = 0;
        Fr = 0;
        Ra = 0;
        Ac = 0;
        Th = 0;
        Pa = 0;
        U  = 0;
        Np = 0;
        Pu = 0;
        Am = 0;
        Cm = 0;
        Bk = 0;
        Cf = 0;
        Es = 0;
        Fm = 0;
        Md = 0;
        No = 0;
        Lr = 0;
        Rf = 0;
        Db = 0;
        Sg = 0;
        Bh = 0;
        Hs = 0;
        Mt = 0;
        rGroups = struct;
    end

    properties (Constant=true)
        elementNames = {
            'H'; 'He'; 'Li'; 'Be';  'B';  'C';  'N';  'O';  'F'; 'Ne'; 'Na'; 'Mg'; 'Al'; 'Si';  'P';  'S';
            'Cl'; 'Ar';  'K'; 'Ca'; 'Sc'; 'Ti';  'V'; 'Cr'; 'Mn'; 'Fe'; 'Co'; 'Ni'; 'Cu'; 'Zn'; 'Ga'; 'Ge';
            'As'; 'Se'; 'Br'; 'Kr'; 'Rb'; 'Sr';  'Y'; 'Zr'; 'Nb'; 'Mo'; 'Tc'; 'Ru'; 'Rh'; 'Pd'; 'Ag'; 'Cd';
            'In'; 'Sn'; 'Sb'; 'Te';  'I'; 'Xe'; 'Cs'; 'Ba'; 'La'; 'Ce'; 'Pr'; 'Nd'; 'Pm'; 'Sm'; 'Eu'; 'Gd';
            'Tb'; 'Dy'; 'Ho'; 'Er'; 'Tm'; 'Yb'; 'Lu'; 'Hf'; 'Ta';  'W'; 'Re'; 'Os'; 'Ir'; 'Pt'; 'Au'; 'Hg';
            'Tl'; 'Pb'; 'Bi'; 'Po'; 'At'; 'Rn'; 'Fr'; 'Ra'; 'Ac'; 'Th'; 'Pa';  'U'; 'Np'; 'Pu'; 'Am'; 'Cm';
            'Bk'; 'Cf'; 'Es'; 'Fm'; 'Md'; 'No'; 'Lr'; 'Rf'; 'Db'; 'Sg'; 'Bh'; 'Hs'; 'Mt'};
    end

    methods
        %4 ways to instatiate:
        %(1) x = EmpiricalFormula()
        %(2) x = EmpiricalFormula(string)
        %(3) x = EmpiricalFormula(cell array of strings)
        %(4) x = EmpiricalFormula(struct)
        function this = EmpiricalFormula(formula)
            import edu.stanford.covert.util.EmpiricalFormula;

            if ~exist('formula','var')
                return;
            end

            if ischar(formula) && any(size(formula))
                formula = {formula};
            end

            this = EmpiricalFormula.empty(numel(formula),0);

            if ischar(formula)
            elseif iscell(formula)
                this(numel(formula),1) = EmpiricalFormula;
                this = reshape(this, size(formula));

                for i=1:numel(formula)
                    tokens = regexp(formula{i},'([A-Z][a-z]*)(\d*)', 'tokens');
                    for j=1:length(tokens)
                        if isempty(tokens{j}{2});
                            tokens{j}{2}=1;
                        else
                            tokens{j}{2}=str2double(tokens{j}{2});
                        end

                        if tokens{j}{2}==0
                            continue;
                        end

                        if ismember(tokens{j}{1}, EmpiricalFormula.elementNames)
                            if this(i).(tokens{j}{1})
                                throw(MException('EmpiricalFormula:invalidEmpiricalFormula','Elements cannot be repeated within a chemical formula'));
                            else
                                this(i).(tokens{j}{1})=tokens{j}{2};
                            end
                        else
                            if isfield(this(i).rGroups, tokens{j}{1}) && this(i).rGroups.(tokens{j}{1})
                                throw(MException('EmpiricalFormula:invalidEmpiricalFormula','Elements cannot be repeated within a chemical formula'));
                            else
                                this(i).rGroups.(tokens{j}{1})=tokens{j}{2};
                            end
                        end
                    end
                end
            elseif isstruct(formula)
                this(numel(formula),1) = EmpiricalFormula;
                this = reshape(this, size(formula));

                elements = fieldnames(formula);
                for i=1:numel(formula)
                    for j=1:numel(elements)
                        if ismember(elements{j}, EmpiricalFormula.elementNames)
                            this(i).(elements{j}) = double(formula(i).(elements{j}));
                        elseif formula(i).(elements{j})
                            this(i).rGroups.(elements{j}) = double(formula(i).(elements{j}));
                        end
                    end
                end
            else
                throw(MException('EmpiricalFormula:unsupportedSyntax','No constructor matches the input'));
            end
        end
    end

    %molecular weight, number
    methods
        function value = weight(this)
            import edu.stanford.covert.util.ConstantUtil;
            import edu.stanford.covert.util.EmpiricalFormula;

            value = zeros(size(this));

            for i=1:numel(this)
                if ~isempty(fieldnames(this(i).rGroups))
                    value(i) = NaN;
                    continue;
                end

                for j=1:numel(EmpiricalFormula.elementNames)
                    value(i) = value(i) + this(i).(EmpiricalFormula.elementNames{j}) * ConstantUtil.elements.(EmpiricalFormula.elementNames{j});
                end
            end
        end

        function value = mass(this)
            value = weight(this);
        end

        function value = number(this)
            import edu.stanford.covert.util.ConstantUtil;
            import edu.stanford.covert.util.EmpiricalFormula;

            value = zeros(size(this));

            for i=1:numel(this)
                if ~isempty(fieldnames(this(i).rGroups))
                    value(i) = NaN;
                    continue;
                end

                for j=1:numel(EmpiricalFormula.elementNames)
                    value(i) = value(i) + this(i).(EmpiricalFormula.elementNames{j}) * ConstantUtil.elementNumbers.(EmpiricalFormula.elementNames{j});
                end
            end
        end
    end

    %operators
    methods
        function tf = isequal(A, B)
            import edu.stanford.covert.util.EmpiricalFormula;

            tf = false;

            if ~isequal(class(A), class(B))
                return;
            end

            if ~isequal(size(A), size(B))
                return;
            end

            for i=1:numel(A)
                for j=1:numel(EmpiricalFormula.elementNames)
                    if ~isequal(A(i).(EmpiricalFormula.elementNames{j}), B(i).(EmpiricalFormula.elementNames{j}))
                        return;
                    end
                end

                rGroupsA = fieldnames(A(i).rGroups);
                rGroupsB = fieldnames(B(i).rGroups);
                rGroupsAB = intersect(rGroupsA,rGroupsB);
                rGroupsA = setdiff(rGroupsA, rGroupsAB);
                rGroupsB = setdiff(rGroupsB, rGroupsAB);
                for j=1:numel(rGroupsAB)
                    if ~isequal(A(i).rGroups.(rGroupsAB{j}), B(i).rGroups.(rGroupsAB{j}))
                        return;
                    end
                end
                for j=1:numel(rGroupsA)
                    if ~isequal(A(i).rGroups.(rGroupsA{j}), 0)
                        return;
                    end
                end
                for j=1:numel(rGroupsB)
                    if ~isequal(B(i).rGroups.(rGroupsB{j}), 0)
                        return;
                    end
                end
            end

            tf = true;
        end

        function C = eq(A,B)
            import edu.stanford.covert.util.EmpiricalFormula;

            if ~isa(A, 'EmpiricalFormula')
                A = EmpiricalFormula(A);
            end

            if ~isa(B, 'EmpiricalFormula')
                B = EmpiricalFormula(B);
            end

            if numel(A)==1
                tmp=A;
                A=B;
                B=tmp;
            end

            C = false(size(A));

            if numel(B)==1
                for i=1:numel(A)
                    C(i) = isequal(A(i), B);
                end
            else
                if ~isequal(size(A), size(B))
                    throw(MException('EmpiricalFormula:invalidDimensions','Sizes of A and B must be equal'));
                end

                for i=1:numel(A)
                    C(i)=isequal(A(i), B(i));
                end
            end
        end

        function C = ne(A,B)
            C = ~eq(A,B);
        end
    end

    %algebra
    methods
        function C = uplus(A)
            C = A;
        end

        function C = uminus(A)
            import edu.stanford.covert.util.EmpiricalFormula;

            C = A;

            for i=1:numel(A)
                for j=1:numel(EmpiricalFormula.elementNames)
                    C(i).(EmpiricalFormula.elementNames{j})=-A(i).(EmpiricalFormula.elementNames{j});
                end

                rGroups = fieldnames(A(i).rGroups);
                for j=1:numel(rGroups)
                    C(i).rGroups.(rGroups{j}) = -A(i).rGroups.(rGroups{j});
                end
            end
        end

        function C = plus(A, B)
            import edu.stanford.covert.util.EmpiricalFormula;

            if ~isa(A, 'EmpiricalFormula')
                A = EmpiricalFormula(A);
            end

            if ~isa(B, 'EmpiricalFormula')
                B = EmpiricalFormula(B);
            end

            if numel(A)==1
                tmp=A;
                A=B;
                B=tmp;
            end

            C = A;

            if numel(B)==1
                rGroups = fieldnames(B.rGroups);

                for i=1:numel(A)
                    for j=1:numel(EmpiricalFormula.elementNames)
                        C(i).(EmpiricalFormula.elementNames{j}) = A(i).(EmpiricalFormula.elementNames{j}) + B.(EmpiricalFormula.elementNames{j});
                    end

                    for j=1:numel(rGroups)
                        if isfield(A(i).rGroups, rGroups{j});
                            C(i).rGroups.(rGroups{j}) = A(i).rGroups.(rGroups{j}) + B.rGroups.(rGroups{j});
                            if ~C(i).rGroups.(rGroups{j})
                                C(i).rGroups = rmfield(C(i).rGroups, rGroups{j});
                            end
                        else
                            C(i).rGroups.(rGroups{j}) = B.rGroups.(rGroups{j});
                        end
                    end
                end
            else
                if ~isequal(size(A), size(B))
                    throw(MException('EmpiricalFormula:invalidDimensions','Sizes of A and B must be equal'));
                end

                for i=1:numel(A)
                    for j=1:numel(EmpiricalFormula.elementNames)
                        C(i).(EmpiricalFormula.elementNames{j}) = A(i).(EmpiricalFormula.elementNames{j}) + B(i).(EmpiricalFormula.elementNames{j});
                    end

                    rGroups = fieldnames(B(i).rGroups);
                    for j=1:numel(rGroups)
                        if isfield(A(i).rGroups, rGroups{j});
                            C(i).rGroups.(rGroups{j}) = A(i).rGroups.(rGroups{j}) + B(i).rGroups.(rGroups{j});
                            if ~C(i).rGroups.(rGroups{j})
                                C(i).rGroups = rmfield(C(i).rGroups, rGroups{j});
                            end
                        else
                            C(i).rGroups.(rGroups{j}) = B(i).rGroups.(rGroups{j});
                        end
                    end
                end
            end
        end

        function C = minus(A, B)
            C = A + -B;
        end

        function C = times(A, B)
            C = A * B;
        end

        function C = mtimes(A, B)
            import edu.stanford.covert.util.EmpiricalFormula;

            if ~isa(A, 'EmpiricalFormula')
                tmp = A;
                A = B;
                B = tmp;
            end

            if ~isnumeric(B) || numel(B)~=1
                throw(MException('EmpiricalFormula:error','EmpiricalFormula only supports scalar multiplication'))
            end

            C = A;

            for i=1:numel(A)
                for j=1:numel(EmpiricalFormula.elementNames)
                    C(i).(EmpiricalFormula.elementNames{j}) = A(i).(EmpiricalFormula.elementNames{j}) * B;
                end

                rGroups = fieldnames(A(i).rGroups);
                for j=1:numel(rGroups)
                    C(i).rGroups.(rGroups{j}) = A(i).rGroups.(rGroups{j}) * B;
                    if ~C(i).rGroups.(rGroups{j})
                        C(i).rGroups = rmfield(C(i).rGroups, rGroups{j});
                    end
                end
            end
        end

        function C = rdivide(A, B)
            C = A / B;
        end

        function C = ldivide(A, B)
            C = B / A;
        end

        function C = mrdivide(A, B)
            import edu.stanford.covert.util.EmpiricalFormula;

            if isa(B, 'EmpiricalFormula')
                throw(MException('EmpiricalFormula:error','Unsupported division'))
            end

            C = A * (1/B);
        end

        function C = mldivide(A, B)
            C = B / A;
        end
    end

    %reduction
    methods
        function C = sum(A)
            import edu.stanford.covert.util.EmpiricalFormula;

            if numel(A)==1
                C=A;
                return;
            end

            C = EmpiricalFormula();

            for i=1:numel(A)
                C = C + A(i);
            end
        end
    end

    %type casting
    methods
        function formula = cast(this, newclass)
            switch newclass
                case 'cell',    formula = cell(this);
                case 'char',    formula = char(this);
                case 'struct',  formula = struct(this);
                otherwise,      throw(MException('EmpiricalFormula:error','unable to convert EmpiricalFormula to %s',newclass));
            end
        end

        function formula = cell(this)
            formula = cell(size(this));
            for i=1:numel(this)
                formula{i}=char(this(i));
            end
        end

        function formula = char(this)
            import edu.stanford.covert.util.EmpiricalFormula;

            if numel(this)>1
                throw(MException('ChemicalFormula:error','only a scalar chemical formula can be converted to a string'));
            end

            formula = '';

            for j=1:numel(EmpiricalFormula.elementNames)
                if this.(EmpiricalFormula.elementNames{j})==1
                    formula = [formula EmpiricalFormula.elementNames{j}];
                elseif this.(EmpiricalFormula.elementNames{j})>1
                    formula = [formula EmpiricalFormula.elementNames{j} num2str(this.(EmpiricalFormula.elementNames{j}))];
                end
            end

            rGroups = fieldnames(this.rGroups);
            for j=1:numel(rGroups)
                if this.rGroups.(rGroups{j})==1
                    formula = [formula rGroups{j}];
                elseif this.rGroups.(rGroups{j})>1
                    formula = [formula rGroups{j} num2str(this.rGroups.(rGroups{j}))];
                end
            end
        end

        function formula = struct(this)
            import edu.stanford.covert.util.EmpiricalFormula;

            formula = repmat(struct, size(this));
            for i=1:numel(this)
                for j=1:numel(EmpiricalFormula.elementNames)
                    if this(i).(EmpiricalFormula.elementNames{j})
                        formula(i).(EmpiricalFormula.elementNames{j})=this(i).(EmpiricalFormula.elementNames{j});
                    end
                end

                rGroups = fieldnames(this(i).rGroups);
                for j=1:numel(rGroups)
                    formula(i).(rGroups{j})=this(i).rGroups.(rGroups{j});
                end
            end
        end
    end

    %display
    methods
        function disp(this)
            this.display();
        end

        function display(this)
            formulae = cell(this);
            fprintf('%dx1 EmpiricalFormula\n',numel(this));
            for i=1:numel(this)
                fprintf('%4d.\t%s\n',i,formulae{i});
            end
        end
    end
end