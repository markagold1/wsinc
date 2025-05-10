function y = inrange(x,lwr,upr,range_type)
% Usage: y = inrange(x,lwr,upr,range_type)
%
% Determine which elements of an array are within the specified range.
%
%   x................numeric array to evaluate
%   lwr..............numeric scalar lower limit of range
%   upr..............numeric scalar upper limit of range
%   range_type.......char array specifying type of range interval
%      'closed'        include upper and lower limits (default)
%      'open'          exclude upper and lower limits
%      'openleft' or    
%      'closedright'   exclude lower and include upper limit
%      'closedleft' or  
%      'openright'     include lower and exclude upper limit
%   y................logical array of results where a "1" in position
%                    n means that the nth element of x falls within
%                    the specified range
%
% Example: Find all elements in x that fall within the closed range [3,5].
%   x = [1,2,3,4,5,6,7,8];
%   y = inrange(x,3,5,'closed')
%   y =
%     0   0   1   1   1   0   0   0
%
    if nargin < 4
        range_type = 'closed';
    end

    switch range_type
        case 'closed'
            y = and(ge(x,lwr), le(x,upr));
        case 'open'
            y = and(gt(x,lwr), lt(x,upr));
        case {'openleft', 'closedright'}
            y = and(gt(x,lwr), le(x,upr));
        case {'openright', 'closedleft'}
            y = and(ge(x,lwr), lt(x,upr));
        otherwise
            y = nan;
            fprintf(2,'Invalid range type.\n');
    end

end % function
