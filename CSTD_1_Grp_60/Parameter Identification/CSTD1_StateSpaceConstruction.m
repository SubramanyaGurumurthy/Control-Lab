
% TUHH :: ICS
% Control Lab CSTD1
% State space construction
% Last update: 10.11.2014
% AMG

function [sys]=CSTD1_StateSpaceConstruction(J1,J2,J3,d1,d2,d3,k1,k2)

%% Functionality: construct the state space model based on identified parameters.


% A matrix Nx-By-Nx: 6x6
A=[  0      1     0           0      0      0     ;
    -k1/J1 -d1/J1 k1/J1       0      0      0     ;
     0      0     0           1      0      0     ;
     k1/J2  0    -(k1+k2)/J2 -d2/J2  k2/J2  0     ;
     0      0     0           0      0      1     ;
     0      0     k2/J3       0     -k2/J3 -d3/J3 ];

% B matrix Nx-by-Nu: 6x1
B=[ 0    ;  
    1/J1 ;
    0    ;  
    0    ;
    0    ;
    0    ];

% C matrix Ny-by-Nx: 3x6
C=[  1     0     0     0     0     0 ;
     0     0     1     0     0     0 ;
     0     0     0     0     1     0 ];
 
% D matrix Ny-by-Nu: 3x1
D=[ 0 ;
    0 ;
    0 ];


%%
sys=ss(A,B,C,D);
