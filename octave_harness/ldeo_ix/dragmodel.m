%======================================================================
%                    D R A G M O D E L . M 
%                    doc: Wed Jan  7 16:23:52 2009
%                    dlm: Wed Jan  7 16:24:07 2009
%                    (c) 2009 A.M. Thurnherr
%                    uE-Info: 11 0 NIL 0 0 72 0 2 4 NIL ofnI
%======================================================================

% CHANGES BY ANT:
%   Jan  7, 2009: - tightened use of exist()

% 
%	dragmodel
%       predict the position of the CTD given the ships position and CTD depth
% 
% DEFINE DRAG MODEL constants
%

if ~exist('p','var')
 p.name='dragmodel';
end

% CTD/LADCP parameter
%  Weight of CTD/LADCP system in kg
p=setdefv(p,'dmodel_package_weight',400);
%  Horizontal area of CTD/LADCP system in m^2
p=setdefv(p,'dmodel_package_area_h',.2);
%  Vertical area of CTD/LADCP system in m^2
p=setdefv(p,'dmodel_package_area_v',0.1);
% drag coefficient for CTD/LADCP system (assume cylinder)
p=setdefv(p,'dmodel_package_cd',1.3);

% Cable parameter
% Cable weight in kg per meter length
p=setdefv(p,'dmodel_wire_weight',0.5);
% Cable wire diameter fraction of inch
p=setdefv(p,'dmodel_wire_diameter_inch',0.5);
% Cable thickness in m is equal to area 
p=setdefv(p,'dmodel_wire_area',0.0254*p.dmodel_wire_diameter_inch);
% Cable drag coefficient
p=setdefv(p,'dmodel_wire_cd',1.3);

% Model parameter:
% length of wire segments
p=setdefv(p,'dmodel_wire_length',10);

% Model design:
%
% each element has an angle with the vertical: alpha 
%    where 0 means pointing straight down
%
% each element had an horizontal angle: beta
%   where 0 means pointing to the north and positive angles point east
% each element has a u,v,w velocity defined at its center
%
% At each junction between elements the forces are calculated
%  in the north/east/up direction units are N=kg*m/s^2
%
% Acceleration of gravity  m/s^2
g=9.81;

% Density of water
rho0=1028;

if ~exist('al','var')
 % initialize model
 ncable=10;
 cable_l=ones(ncable,1)*p.dmodel_wire_length;
 % angle to the vertical
 al=zeros(ncable+1,1);
 % angle in horizontal
 be=zeros(ncable+1,1);
 % node velocity
 u=zeros(ncable+1,1);
 v=zeros(ncable+1,1);
 % node acceleration
 dudt=zeros(ncable+1,1);
 dvdt=zeros(ncable+1,1);
 % node position
 X=zeros(ncable+1,1);
 Y=zeros(ncable+1,1);
 Z=zeros(ncable+1,1);
 % node horizontal distance
 H=zeros(ncable+1,1);
 % node vertical distance
 D=zeros(ncable+1,1);
end

if ~exist('dr','var')
 dr.z=[0 3000];
 dr.u=[0 0];
 dr.v=[0 0];
 d.wm=[0 0];
end


% time step
nt=1;
dt=60;

% precalculate 
sal=sin(al*pi/180);
cal=cos(al*pi/180);
sbe=sin(be*pi/180);
cbe=cos(be*pi/180);

% get depth of each node
Z=al*0;
for i=1:length(cable_l)
 in=length(cable_l)-i+1;
 Z(in)=Z(in+1)-p.dmodel_wire_length*cos(al(in));
end

% velocity of nodes
vel=u+sqrt(-1)*v;

% interpolate ocean velocity to drag model points
uo=interp1(-dr.z,dr.u,Z-p.dmodel_wire_length/2,'linear','extrap');
vo=interp1(-dr.z,dr.v,Z-p.dmodel_wire_length/2,'linear','extrap');
velo=uo+sqrt(-1)*vo;
wctd=d.wm(nt);

%  Compute forces from bottom of CTD up to the ship
% 
% CTD/LADCP first
% gravity
Fv(1)=-p.dmodel_package_weight*g;

% drag
% normal drag Force in direction of flow
Fdn(1)=(-sal(1)*wctd+cal(1)*abs(vel(1)-velo(1))) * ...
       abs(sal(1)*wctd+cal(1)*abs(vel(1)-velo(1))) * ...
      p.dmodel_package_area_h*p.dmodel_package_cd*rho0;
% along drag Force
Fda(1)=(-cal(1)*wctd+sal(1)*abs(vel(1)-velo(1))) * ...
       abs(cal(1)*wctd+sal(1)*abs(vel(1)-velo(1))) * ...
      p.dmodel_package_area_h*p.dmodel_package_cd*rho0;

% project in earth coordinates
% vertical drag
Fv(1)=Fv(1)+Fdn(1)*sal(1)+Fda(1)*cal(1);

% tension
Ft(1)=Fda(1)*sal(1)+Fv(1)*cal(1);

% horizontal drag
Fh(1)=Fdn(1)*cal(1)+Ft(1)*sal(1);
Fhx(1)=Fh(1)*cos(angle(vel(1)-velo(1)));
Fhy(1)=Fh(1)*sin(angle(vel(1)-velo(1)));

% loop over the wire elements to calculate Forces and Tension
for i=1:length(cable_l)
 i1=i+1;
 % gravity
 Fv(i1)=Fv(i)-p.dmodel_wire_weight*p.dmodel_wire_length*g;

 % drag
 % normal drag Force in direction of flow
 Fdn(i1)=(-sal(i1)*wctd+cal(i1)*abs(vel(i1)-velo(i1))) * ...
         abs(sal(i1)*wctd+cal(i1)*abs(vel(i1)-velo(i1))) * ...
      p.dmodel_wire_area*p.dmodel_wire_cd*p.dmodel_wire_length;
 % along drag Force
 Fda(i1)=(-cal(i1)*wctd+sal(i1)*abs(vel(i1)-velo(i1))) * ...
         abs(cal(i1)*wctd+sal(i1)*abs(vel(i1)-velo(i1))) * ...
      p.dmodel_wire_area*p.dmodel_wire_cd*p.dmodel_wire_length;

 % project in earth coordinates
 % vertical drag
 Fv(i1)=Fv(i1)+Fdn(i1)*sal(i1)+Fda(i1)*cal(i1);

 % tension
 Ft(i1)=Fda(i1)*sal(i1)+Fv(i1)*cal(i1);

 % horizontal drag
 Fh(i1)=Fdn(i1)*cal(i1)+Ft(i1)*sal(i1);
 Fhx(i1)=Fh(i1)*cos(angle(vel(i1)-velo(i1)));
 Fhy(i1)=Fh(i1)*sin(angle(vel(i1)-velo(i1)));

end

% loop over wire elements to calculate accelation from F = m*a
%  a = F/m
% use "dynamical mass" Fv/g

dudt(1) = Fhx(1)/abs(Fv(1)/g);
dvdt(1) = Fhy(1)/abs(Fv(1)/g);


for i=1:(length(cable_l))
  i1=i+1;
  dudt(i1) = 0.5*(Fhx(i1)+Fhx(i))/abs(Fv(i1)/g);
  dvdt(i1) = 0.5*(Fhy(i1)+Fhy(i))/abs(Fv(i1)/g);
end 

% loop to update velocity of nodes

for i=0:(length(cable_l))
  i1=i+1;
  u(i1) = u(i1) + dudt(i1)*dt/1000;
  v(i1) = v(i1) + dvdt(i1)*dt/1000;
end 


% loop to calculate horizontal position of nodes
for i=0:(length(cable_l))
  i1=i+1;
  X(i1)=X(i1)+u(i1)*dt;
  Y(i1)=Y(i1)+v(i1)*dt;
end

% loop to calculate vertical position of nodes
Z(length(cable_l)+1)=0;
for i=1:(length(cable_l))
  in=length(cable_l)+1-i;
  in1=in+1;
  H(in1)=sqrt((X(in1)-X(in)).^2 + (Y(in1)-Y(in)).^2);
  D(in1)=sqrt(p.dmodel_wire_length.^2 -  H(in).^2);
  Z(in)=Z(in1)-D(in1);
  al(in)=atan2(H(in),D(in))*180/pi;
  be(in)=atan2(Y(in1)-Y(in),X(in1)-X(in))*180/pi;
end

