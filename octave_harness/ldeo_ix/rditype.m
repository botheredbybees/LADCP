%======================================================================
%                    R D I T Y P E . M 
%                    doc: Wed Jan  7 16:49:09 2009
%                    dlm: Wed Jan  7 16:49:22 2009
%                    (c) 2009 A.M. Thurnherr
%                    uE-Info: 287 0 NIL 0 0 72 0 2 8 NIL ofnI
%======================================================================

% CHANGES BY ANT:
%   Jan  7, 2009: - tightened use of exist()

function [data,d] = rditype(file, verbose, select);

% function data = rditype(file, verbose, select);
%	Read the fixed leader data from a raw ADCP
%	data file.
%	Returns the contents of the fixed leader
%	as 31 elements of the vector 'data' or an
%	empty matrix if the fixed leader ID is not
%	identified (error condition)
%	If the variable select is provided as a vector
%	of zeros and ones, the function will return
%	only the elements of data which correspond to
%	a one in the vector select.  Select must be the
%	the same length as the number of fields in the
%	record, currently 32.
%	Set verbose=1 for a text output.

% get fixed leader data
 if nargin<3, select = ones(1,31); end
 if nargin<2, verbose = 0; end
 if nargin<1, disp(' provide file name '), return, end

 
 fid = fopen(file,'r','l');
 hid = 127;  % header identification byte
 [id,n] = fread(fid,2,'uint8');
 status = fseek(fid,0,'bof');
 if id(1) ~= hid
  % Fix parameter for NB-ADCP Kiel
  data=zeros(1,32);
  d.nbytes=NaN;
  d.ntypes=NaN;
  d.Firm_Version=0;
  d.Frequency=150;
  d.Up=0;
  d.Beam_angle=20;
  d.Cell_length=1600;
  d.Pings_per_Ensemble=8;
  d.Pulse_length=1600;
  d.Blank=NaN;
  d.Mode=NaN;
  d.Time_Pings=NaN;
  d.Coordinates=3;
  d.NarrowBand=1;
  return
 end
 [d.nbytes, d.ntypes, offsets] = rdhead(fid, verbose);

 % ---- now at start of fixed leader data

% originally 
% Written by Marinna Martini
% for the U.S. Geological Survey
% Atlantic Marine Geology, Woods Hole, MA
% 1/7/95

data=zeros(1,32);
fld=1;  
% make sure we're looking at the beginning of
% the fixed leader record by testing for it's ID
data(fld)=fread(fid,1,'ushort');
if(data(fld)~=0),
	disp('Fixed Leader ID not found');
	data=[];
	return;
end
fld=fld+1;
% version number of CPU firmware
data(fld)=fread(fid,1,'uchar');
fld=fld+1;
% revision number of CPU firmware
data(fld)=fread(fid,1,'uchar');
d.Firm_Version=data(fld-1)+data(fld)/100;
if verbose, disp(sprintf('CPU Version %d.%d',data(fld-1),data(fld))); end;
fld=fld+1;
% configuration, uninterpreted
data(fld)=fread(fid,1,'uchar');
	b=dec2binv(data(fld));
	freqs=[75 150 300 600 1200 2400];
	junk=bin2decv(b(6:8));
d.Frequency=freqs(junk+1);
d.Up=str2num(b(1));
if verbose, 
	disp(sprintf('Hardware Configuration for LSB %d',data(fld))); 
	disp(sprintf('	System Frequency = %d kHz',freqs(junk+1))); 
	if b(5) == '0', disp('	Concave Beam'); end		
	if b(5) == '1', disp('	Convex Beam'); end		
	junk=bin2decv(b(3:4));
	disp(sprintf('Sensor Configuration #%d',junk+1)); 
	if b(2) == '0', disp('	Transducer head not attached'); end		
	if b(2) == '1', disp('	Transducer head attached'); end		
	if b(1) == '0', disp('	Downward facing beam orientation'); end		
	if b(1) == '1', disp('	Upward facing beam orientation'); end		
end;
fld=fld+1;
data(fld)=fread(fid,1,'uchar');
	b=dec2binv(data(fld));
	angles = [15 20 30 0];
	junk=bin2decv(b(7:8));
d.Beam_angle=angles(junk+1);
if verbose, 
	disp(sprintf('Hardware Configuration MSB %d',data(fld))); 
	disp(sprintf('	Beam angle = %d degrees',angles(junk+1))); 
	junk=bin2decv(b(1:4));
	if junk == 4, disp('	4-beam janus configuration'); end
	if junk == 5, disp('	5-beam janus configuration, 3 demodulators'); end
	if junk == 15, disp('	4-beam janus configuration, 2 demodulators'); end
end;
fld=fld+1;
% real (0) or simulated (1) data flag
data(fld)=fread(fid,1,'uchar');	fld=fld+1;
% undefined
data(fld)=fread(fid,1,'uchar');	fld=fld+1;
% number of beams
data(fld)=fread(fid,1,'uchar');	fld=fld+1;
% number of depth cells
data(fld)=fread(fid,1,'uchar');
d.Depth_Cells=data(fld);
if verbose, disp(sprintf('Number of depth cells %d',data(fld))); end;
fld=fld+1;
% pings per ensemble
data(fld)=fread(fid,1,'ushort');
d.Pings_per_Ensemble=data(fld);
if verbose, disp(sprintf('Pings per ensemble %d',data(fld))); end;
fld=fld+1;
% depth cell length in cm
data(fld)=fread(fid,1,'ushort');
d.Cell_length=data(fld);
if verbose, disp(sprintf('Depth cell size %d cm',data(fld))); end
fld=fld+1;
% blanking distance (WF)
data(fld)=fread(fid,1,'ushort');
d.Blank=data(fld);
if verbose, disp(sprintf('Blank after xmit distance %d cm',data(fld))); end
fld=fld+1;
% Profiling mode (WM)
data(fld)=fread(fid,1,'uchar');
d.Mode=data(fld);
if verbose, disp(sprintf('Profiling mode %d',data(fld))); end
fld=fld+1;
% Minimum correlation threshold (WC)
data(fld)=fread(fid,1,'uchar');
d.Min_Correlation=data(fld);
if verbose, disp(sprintf('Correlation threshold %d',data(fld))); end
fld=fld+1;
% number of code repetitions
data(fld)=fread(fid,1,'uchar');
d.Code_rep=data(fld);
fld=fld+1;
% Minimum percent good to output data (WG)
data(fld)=fread(fid,1,'uchar');
d.Min_Percgood=data(fld);
fld=fld+1;
% Error velocity threshold (WE)
d.Max_Errorvel=data(fld);
data(fld)=fread(fid,1,'ushort');
if verbose, disp(sprintf('Error Velocity Threshold %d mm/s',data(fld))); end
fld=fld+1;
% time between ping groups (TP)
data(fld)=fread(fid,1,'uchar');
d.Time_Pings=data(fld)*60;
fld=fld+1;
data(fld)=fread(fid,1,'uchar');
d.Time_Pings=d.Time_Pings+data(fld);
fld=fld+1;
data(fld)=fread(fid,1,'uchar');
d.Time_Pings=d.Time_Pings+data(fld)/100;
if verbose, disp(sprintf('Time between ping groups %d:%d.%d',...
	data(fld-2), data(fld-1), data(fld))); end
fld=fld+1;
% coordinate transformation (EX)
data(fld)=fread(fid,1,'uchar');
	b=dec2binv(data(fld));
	junk=bin2decv(b(4:5));
d.Coordinates=junk;
d.use_tilt=str2num(b(6));
if verbose, 
	disp(sprintf('Coordinate Transformation = %d',data(fld))); 
	if junk == 0, disp('	Data stored coordinates = Beam'); end
	if junk == 1, disp('	Data stored coordinates = Instrument'); end
	if junk == 2, disp('	Data stored coordinates = Ship'); end
	if junk == 3, disp('	Data stored coordinates = Earth'); end
	if b(6) == '1', disp('	Tilts used in transformation'); end
%	if b(7) == '1', disp('	3-beam solution used, this ensemble'); end
end
fld=fld+1;
% Heading Alignment (EA)
data(fld)=fread(fid,1,'uint16');
d.headin_alignment=data(fld);
fld=fld+1;
% Heading Bias (EB)
data(fld)=fread(fid,1,'uint16');
d.headin_bias=data(fld);
if verbose, disp(sprintf('Heading Bias: %d deg',data(fld)./100)); end
fld=fld+1;
% Sensor source (EZ)
data(fld)=fread(fid,1,'uchar');
	b=dec2binv(data(fld));
d.sensor_source=b;
if verbose,
	disp(sprintf('Sensor Source = %d',data(fld))); 
	if b(2) == '1', disp('	Sound speed computed from ED, ES, ET'); end
	if b(3) == '1', disp('	ED taken from depth sensor'); end	
	if b(4) == '1', disp('	EH taken from xducer heading sensor'); end	
	if b(5) == '1', disp('	EP taken from xducer pitch sensor'); end	
	if b(6) == '1', disp('	ER taken from xducer roll sensor'); end	
	if b(7) == '1', disp('	ES derived from conductivity sensor'); end	
	if b(8) == '1', disp('	ET taken from temperature sensor'); end	
end
fld=fld+1;
% Sensors available
data(fld)=fread(fid,1,'uchar');
	b=dec2binv(data(fld));
d.sensor_avail=b;
if verbose,
	disp(sprintf('Sensor Availability = %d',data(fld))); 
	if b(3) == '1', disp('	depth sensor'); end	
	if b(4) == '1', disp('	heading sensor'); end	
	if b(5) == '1', disp('	pitch sensor'); end	
	if b(6) == '1', disp('	roll sensor'); end	
	if b(7) == '1', disp('	conductivity sensor'); end	
	if b(8) == '1', disp('	temperature sensor'); end	
end
fld=fld+1;
% Bin 1 distance
data(fld)=fread(fid,1,'ushort');
if verbose, disp(sprintf('Distance to center of bin 1: %d cm',data(fld))); end
fld=fld+1;
% xmit pulse length
data(fld)=fread(fid,1,'ushort');
d.Pulse_length=data(fld);
fld=fld+1;
% starting depth cell
data(fld)=fread(fid,1,'uchar');
fld=fld+1;
% ending depth cell
data(fld)=fread(fid,1,'uchar');
fld=fld+1;
% false target reject threshold
data(fld)=fread(fid,1,'uchar');
d.target_max=data(fld);
fld=fld+1;
% spare
data(fld)=fread(fid,1,'uchar');
fld=fld+1;
% transmit lag distance
data(fld)=fread(fid,1,'ushort');
d.xmit_lag=data(fld);

if length(select) == length(data),
	data(find(select==0))=[];
end


%--------------------------------------


function [nb, nt, off] = rdhead(fid, verbose);

%function [nb, nt, off] = rdhead(fid, verbose);
%	Read the header data from a raw ADCP
%	data file opened for binary reading.
%	fid = file handle returned by fopen
%	nb = number of bytes in the ensemble
%	nt = number of data types
%	off = offset to the data for each type
%	Set verbose = 1 for a text output.

% Written by Marinna Martini
% for the U.S. Geological Survey
% Atlantic Marine Geology, Woods Hole, MA
% 1/7/95

data = zeros(1,2);
fld=1;
if ~exist('verbose','var'),
	verbose = 0;
end
nb=[]; nt=[]; off=[];

% make sure we're looking at the beginning of
% the header record by testing for it's ID
junk=fread(fid,2,'uint8');
if((length(junk)~=2) | (ftell(fid)<0)),
	disp('End of file found in rdhead.');
	return;
end
if ((junk(1)~=127) | (junk(2)~=127)),
	disp('Header ID not found');
	return;
end
% get the number of bytes this ensemble
nb = fread(fid,1,'uint16');
if verbose, disp(sprintf('Number of bytes per ensemble %d',nb)); end;
% get the number of data types
fseek(fid,1,'cof');	% skip spare byte position
nt=fread(fid,1,'uchar');
if verbose, disp(sprintf('Number of data types %d',nt)); end;
% get the type offset
off=zeros(nt,1);
for j=1:nt, off(j)=fread(fid,1,'uint16'); end

% ---------------------------------

function d = bin2decv(h)
%BIN2DEC  BIN2DEC('X') returns the decimal number corresponding to the
%        binary number in quotes.  
%		 For example, BIN2DEC([1 1 0 0]) returns 12.

n=length(h);
h=fliplr(h);
d=0;
for i=1:n
	if isstr(h(i)),
		p(i) = str2num(h(i)).*(2^(i-1));
	else
		p(i) = h(i).*(2^(i-1));
	end
end
d=sum(p);

% ---------------------------------

function h = dec2binv(d)
%DEC2BIN DEC2BIN(d) returns the binary number corresponding to the decimal
%        number d.  For example, DEC2BIN(202) returns '11001010'.
%
%	
h=dec2bin(d);
bits=length(h)/8;
if bits ~= fix(bits),
	% not an even multiple of 8, so pad
	npad = 8-((bits-fix(bits))*8);
    h0(1:npad)='0';
	h = [h0,h];
end
	
	
