function [ rotationalMatrix ] = getRotationMatrix( axis, angle )
%getRotationMatrix Returns the rotation matrix of "angle" around "axis".
%
%   [ rotationalMatrix ] = getRotationMatrix( axis, angle )
%
%   "angle" must be in radiants while the possible axes are 'x', 'y', 'z'.

S = sin(angle);
C = cos(angle);

switch axis
    case 'x'
        R = [1  0 0 ;
             0  C S ;
             0 -S C];
    case 'y'
        R = [C 0 -S ;
             0 1  0 ;
             S 0  C];
    case 'z'
        R = [ C S 0 ;
             -S C 0 ;
              0 0 1];
end

rotationalMatrix = R;

end

