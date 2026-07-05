function cc=sounds(P0,T,S)
%function cc=sounds(P0,T,S)
% This function is the translation of the FORTRAN REAL FUNCTION SVEL from
% page 49 from "Unesco technical papers in marine science 44"
%
% It returns the sound speed in seawater in meters per second
% Reference: Chen and Millero 1977, JASA, 62, 1129-1135
% Units:
%      Salinity        S   (PSS-78)
%      Temperature     T   Degrees Celsius (IPTS-68)
%      Pressure        P0  Decibars
%
% Returns:
%      Sound Speed     Meters / Second
% Checkvalues:
%  SVEL=1731.995 :Salinity=40.0, Temp.=40.0, Pres.=10000.0

%--- Scale Pressure to bars ---.*/
    P=P0/10.0;
    SR=sqrt(abs(S));
%--- S.*.*2 Term ---.*/
    D=1.727e-3-7.8936e-6.*P;
%--- S.*.*3/2 Term ---.*/
    B1=7.3637e-5+1.7945e-7.*T;
    B0=-1.922e-2-4.42e-5.*T;
    B=B0+B1.*P;
%--- S.*.*1 Term ---.*/
    A3=(-3.389e-13.*T+6.649e-12).*T+1.100e-10;
    A2=((7.988e-12.*T-1.6002e-10).*T+9.1041e-9).*T-3.9064e-7;
    A1=(((-2.0122e-10.*T+1.0507e-8).*T-6.4885e-8).*T-1.2580e-5).*T+9.4742e-5;
    A0=(((-3.21e-8.*T+2.006e-6).*T+7.164e-5).*T-1.262e-2).*T+1.389;
    A=((A3.*P+A2).*P+A1).*P+A0;
%--- S.*.*0 Term ---.*/
    C3=(-2.3643e-12.*T+3.8504e-10).*T-9.7729e-9;
    C2=(((1.0405e-12.*T-2.5335e-10).*T+2.5974e-8).*T-1.7107e-6).*T+3.1260e-5;
    C1=(((-6.1185e-10.*T+1.3621e-7).*T-8.1788e-6).*T+6.8982e-4).*T+0.153563;
    C0=((((3.1464e-9.*T-1.47800e-6).*T+3.3420e-4).*T-5.80852e-2).*T+5.03711).*T+1402.388;
    C=((C3.*P+C2).*P+C1).*P+C0;

    cc=(C+(A+B.*SR+D.*S).*S);
