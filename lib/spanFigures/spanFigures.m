function spanFigures( varargin )
% 
% Spans the figueres on the entire screen
%
% spanFigures() will arrange (tile) all matlab figures 
% spanFigures( hFig ) arranges all figuers in the handle vector hFig

    if nargin == 0
        hFig = get( 0,'children' );
    else
        hFig = varargin{1};
    end
    
    N = numel(hFig);
    
    C = ceil( sqrt(N) );
    R = ceil( N/C );
    set( hFig, 'units', 'normalized' );
    % set( hFig, 'toolbar', 'none' );    
    % set( hFig, 'toolbar', 'figure' );
    w = 1/C;
    h = 0.99/R;
    
    for n = 1 : N
        
        [r c] = ind2sub( [R C], n );
        p = [ (c-1)*w 0.01+(r-1)*h,  w*1, h*0.82 ];
        set( hFig(n),'position',p );     
        figure( hFig(n) );
        
    end