classdef mock < handle & dynamicprops
    %MOCK A handle class (can be passed by reference) that behaves like a struct.

    methods
        function this = mock(varargin)
            o = struct(varargin{:});
            names = fieldnames(o);
            for i = 1:length(names);
                name = names{i};
                addprop(this, name);
                this.(name) = this.convert(o.(name));
            end
        end
        
        function this = subsasgn(this, s, value)
            n = length(s);
            if s(1).type(1) == '.'
                name = s(1).subs;
                if isempty(findprop(this, name))
                    addprop(this, name);
                    if n > 1
                        this.(name) = struct;
                    end
                end
            end
            this = builtin('subsasgn', this, s, this.convert(value));
        end
    end

    methods (Access = private)
        function val = convert(this, val)
            if isa(val, 'function_handle')
                val = this.bind(val);
            end
        end
        
        function boundFn = bind(this, fn)
            function result = bound(varargin)
                result = fn(this, varargin{:});
            end
            boundFn = @bound;
        end
    end
    
end
