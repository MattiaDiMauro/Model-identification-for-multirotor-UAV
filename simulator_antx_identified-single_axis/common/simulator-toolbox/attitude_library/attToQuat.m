function quat = attToQuat( A )
%   Converts the attitude matrix A to a quaternion vector using a modification 
%   of the "Shepperd's Algorithm"
%
%   quat = attToQuat( A )
%
%   The attitude matrix A represents the rotational matrix from the inertial
%   reference frame to the body frame of the vehicle.
%   The quaternion (quat) will be returned as a 4-by-1 vector with the scalar
%   element as the fourth component.
%
%   References:
%	[1] Markley, F. Landis. "Unit quaternion from rotation matrix." 
%       Journal of guidance, control, and dynamics 31.2 (2008): 440.

vector = [A(1,1), A(2,2), A(3,3), trace(A)];
[M, i] = max(vector);

switch i
    case 1
        q = [ 1 + A(1,1) - A(2,2) - A(3,3);
            A(1,2) + A(2,1);
            A(1,3) + A(3,1);
            A(2,3) - A(3,2)];
        
    case 2
        q = [A(2,1) + A(1,2);
            1 + A(2,2) - A(3,3) - A(1,1);
            A(2,3) + A(3,2);
            A(3,1) - A(1,3)];
        
    case 3
        q = [A(3,1) + A(1,3);
            A(3,2) + A(2,3);
            1 + A(3,3) - A(1,1) - A(2,2);
            A(1,2) - A(2,1)];
        
    case 4
        q = [A(2,3) - A(3,2);
            A(3,1) - A(1,3);
            A(1,2) - A(2,1);
            1 + A(1,1) + A(2,2) + A(3,3);];
        
    otherwise
        error('Check the Attitude matrix');
end

quat = q / norm(q);
