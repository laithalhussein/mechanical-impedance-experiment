function h=display_xy_error(x, y, xe, ye, varargin)
for i=1:length(x)
    
    if ~isempty(xe)
        h=plot( [x(i)-xe(i), x(i)+xe(i)], [y(i), y(i)], varargin{:});
        plot( [x(i),x(i)], [y(i)-ye(i), y(i)+ye(i)], varargin{:});
    else
        h=plot( [x(i),x(i)], [y(i)-ye(i), y(i)+ye(i)], varargin{:});
    end
end

end