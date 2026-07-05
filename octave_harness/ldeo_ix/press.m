function p = press(z)
%function p = press(z)
%C-------------------------------------------------------------------------------
%C
%C                            Function = PRESS
%C
%C     This computes pressure in decibars, given the depth in meters.           *
%C     Formula is from the GEOSECS operation group (See computer routine
%C     listings in the El Nino Watch Data Reports).
%C
%C-------------------------------------------------------------------------------

C1=  2.398599584e05;
C2=  5.753279964e10;
C3=  4.833657881e05;
ARG = C2 -  C3 * z;
p = C1 - sqrt( ARG );

