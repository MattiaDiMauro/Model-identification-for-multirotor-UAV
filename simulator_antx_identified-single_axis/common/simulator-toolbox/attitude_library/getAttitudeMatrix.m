function [ A ] = getAttitudeMatrix( axis, angle )
%	Returns the attitude matrix of "angle" around "axis".
%
%   [ A ] = getAttitudeMatrix( axis, angle )
%
%   "angle" must be in radiants while the possible axes are 'x', 'y', 'z'.
%
%   References:
%   [1] Dreier, Mark E. Introduction to helicopter and tiltrotor simulation. 
%       AIAA (American Institute of Aeronautics & Ast, 2007.

S = sin(angle);
C = cos(angle);

switch axis
    case 'x'
        A = [1  0 0 ;
             0  C S ;
             0 -S C];
    case 'y'
        A = [C 0 -S ;
             0 1  0 ;
             S 0  C];
    case 'z'
        A = [ C S 0 ;
             -S C 0 ;
              0 0 1];
end

end

