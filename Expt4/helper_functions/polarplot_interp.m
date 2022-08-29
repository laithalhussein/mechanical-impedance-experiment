function varargout = polarplot_interp(varargin)

theta_spacing = 0.01;    
if nargin>=2
        theta = varargin{1}(:); theta1 = [];
        r = varargin{2}(:);     r1=[];   
        for k=1:length(theta)-1,
            theta_interp = [[theta(k):.01:theta(k+1)]'; theta(k+1)];
            r_interp = interp1(theta(k:k+1),r(k:k+1),theta_interp);
            theta1 = [theta1; theta_interp];
            r1 = [r1; r_interp];
        end
    [varargout{1:nargout}] = polarplot(theta1, r1, varargin{3:end});
end

end