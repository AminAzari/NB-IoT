
%Traffic
N= 10000; % nodes
S=12; %1 packet per 2 hours
p=0.8; % percentage of uplink requests

l1=1000; % bits in sc1&2
m1=10000; % bits in sc1%600 bits in sc1



u=0.002; % sec
tau=0.005; % sec
lambdab=5; % per sec
b=0.2; % fraction of time in which reference signals are scheduled
bt=0.01;

%Power
Pt=0.2; %watt
Pc=0.01; %watt
PI=0.01; %watt
Pl=0.1; % Watt
al=1;

%Coverage
C=2; %number of classes
c=1; %number of repetitions
c2=2; %number of repetitions
MCR=5; % uplink rate in Kbps
MBR=15; % downlink rate in Kbps
ut=l1/MCR; %in msec
dt=floor(m1/MBR); %in msec
TTI=0.001;
f1=0.5;f2=0.5;
 
%Other
M=16; 
Tth=2;
E0=1000; % joule
Dsy1=0.33; % sec
Dsy2=0.66; % sec