function rms = myRMS(x, varargin)
rms=sqrt(nanmean(x.^2, varargin{:}));
return