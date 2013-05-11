% VEC       = rude(LEN,VAL)
%		run-length DEcoding
%
% [LEN,VAL] = rude(VEC)
%		run-length ENcoding
%
% P         = rude;
%		retrieve subroutine handles in structure P
%   P.d		for run-length DEcoding
%			VEC = P.d(LEN,VAL)
%   P.e		for run-length ENcoding
%			[LEN,VAL] = P.e(VEC)
%
% LEN	: repeat each VAL corresponding LEN times to create VEC
% VAL	: 1xN template array of VALs to be repeated LEN times
%	  - numericals
%	  - strings
%	  - cells (any contents)
% VEC	: DEcode = reconstruced output vector from LEN/VAL
%	  ENcode = input vector to be encoded into LEN/VAL
%
% NOTE
% 1:	LEN <= 0 will remove corresponding VALs
% 2:	<NaN>s and <+Inf>s are treated as not-equal (expected behavior!)
%
% USAGE EXAMPLE
%	vec = rude([1 2 3 0 4],[10 inf nan pi 20])
% %	vec = 10 Inf Inf NaN NaN NaN 20 20 20 20
%	[len,val] = rude(vec)
% %	len = 1 1 1 1 1 1 4 % note nan~=nan / inf~=inf!
% %	val = 10 Inf Inf NaN NaN NaN 20
%
%	s = rude;
%	v.x = pi;
%	w.x = pi; % note <v> and <w> are equal!
%	vec = s.d([1 0 3 2 2 2 3],{'a' 'b' 'cd' 1:3 v w magic(3)})
% %	vec = 'a' 'cd' 'cd' 'cd' 1x3D 1x3D 1x1S 1x1S 1x1S 1x1S 3x3D
%	[len,val] = s.e(vec)
% %	len = 1 3 2 4 3
% %	val = 'a' 'cd' 1x3D 1x1S 3x3D

% created:
%	us	18-Nov-2004
% modified:
%	us	30-Nov-2004 20:42:03	/ TMW FEX

function	[p1,p2]=rude(varargin)

		p2=[];
	if	~nargin & nargout
		p1.d=@rl_decode;
		p1.e=@rl_encode;
	elseif	~nargin
		help(mfilename);
		return;
	else
	if	nargin == 1
		[p1,p2]=rl_encode(varargin{1});
	elseif	nargin >= 2
		p1=rl_decode(varargin{1:2});
	end
	end
		return;
%--------------------------------------------------------------------------------
% run-length decoder
function	vec=rl_decode(len,val)

		lx=len>0 & ~(len==inf);
	if	~any(lx)
		vec=[];
		return;
	end
	if	numel(len) ~= numel(val)
		error(...
		sprintf(['rl-decoder: length mismatch\n',...
			 'len = %-1d\n',...
			 'val = %-1d'],...
			  numel(len),numel(val)));
	end
		len=len(lx);
		val=val(lx);
		val=val(:).';
		len=len(:);
		lc=cumsum(len);
		lx=zeros(1,lc(end));
		lx([1;lc(1:end-1)+1])=1;
		lc=cumsum(lx);
		vec=val(lc);
		return;
%--------------------------------------------------------------------------------
% run-length encoder
function	[len,val]=rl_encode(vec)

	switch	class(vec)
	case	'cell'
		[len,val]=rl_encodec(vec);
	case	'char'
		[len,val]=rl_encoden(double(vec));
		val=char(val);
	otherwise
		[len,val]=rl_encoden(vec);
	end
		return;
%--------------------------------------------------------------------------------
% run-length encode doubles
function	[len,val]=rl_encoden(vec)

	if	isempty(vec)
		len=0;
		val=[];
		return;
	end
		vec=vec(:).';
		vx=[1 diff(double(vec))];
		vx=vx~=0;
		val=vec(vx);
		vc=cumsum(vx).';
		len=accumarray(vc,ones(size(vc))).';
		return;
%--------------------------------------------------------------------------------
% run-length encode cells
function	[len,val]=rl_encodec(vec)

		len=0;
		val={};
		vl=length(vec)+1;
		cl=cellfun('length',vec);
		tix=0;
		tlen=1;
	for	i=1:length(cl)
	if	cl(i)
		tmpl=vec{i};
	for	j=i+1:vl
	if	j>length(cl) || ~isequalwithequalnans(tmpl,vec{j})
		tix=tix+1;
		val{tix}=tmpl;
		len(tix)=tlen;
		cl(i:j-1)=0;
		tlen=1;
		break;
	else
		tlen=tlen+1;
	end
	end
	end
	end
		return;
%--------------------------------------------------------------------------------
