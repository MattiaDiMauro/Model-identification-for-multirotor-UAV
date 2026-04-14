%% Tiltrotor mixer
clear all
close all
clc

syms b sigma real

%[F_z L M N]' = map_... * [f_1_v f_2_v f_3_v f_4_v]'

%[f_1_v f_2_v f_3_v f_4_v]' = mixer_... *[F_z L M N]'

% plus configuration (+)
map_plus = [  -1    -1     -1    -1   ;
               0    -b      0     b   ;
               b     0     -b     0   ;
            -sigma sigma -sigma sigma];
 
mixer_plus = inv(map_plus);
mixer_plus = simplify(mixer_plus);

% cross configuration (X)
map_cross = [    -1            -1           -1           -1      ;
             sqrt(2)/2*b  -sqrt(2)/2*b -sqrt(2)/2*b  sqrt(2)/2*b ;
             sqrt(2)/2*b   sqrt(2)/2*b -sqrt(2)/2*b -sqrt(2)/2*b ;
               -sigma         sigma       -sigma        sigma   ];

mixer_cross = inv(map_cross);
mixer_cross = simplify(mixer_cross);

%% From physical to PX4
map_plus_px4 = [map_plus(:,2), map_plus(:,4), map_plus(:,1), map_plus(:,3)];
map_cross_px4 = [map_cross(:,2), map_cross(:,4), map_cross(:,1), map_cross(:,3)];

mixer_plus_px4 = [mixer_plus(2,:) ;
                  mixer_plus(4,:) ;
                  mixer_plus(1,:) ;
                  mixer_plus(3,:)];

mixer_cross_px4 = [mixer_cross(2,:) ;
                   mixer_cross(4,:) ;
                   mixer_cross(1,:) ;
                   mixer_cross(3,:)];
