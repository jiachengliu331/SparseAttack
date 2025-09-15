function [f1,error,objective_value] = pflow(x)

%global VD PGS emission ploss qloss Lind_worst VI fuelcost fuelvlvcost
%x = [44.9644	 18.3171	 10.0001	 10.0003	 12.0001	 1.1000	 1.0461	 1.0570	 1.0606	 1.0516	 1.1000	 3.6165	 0.0590	 4.0508	 4.9874	 4.5766	 4.9971	 3.1509	 4.9787	 2.2571	 1.0423	 0.9008	 1.0245	 0.9790];
%x = [48.7 21.45 21.06 11.96 12 1.08 1.06 1.03 1.03 1.09 1.04 2.37 2.58 4.2 5 3.7 4.96 3.08 4.99 2.49 1.05 0.94 0.97 0.975];
%x = [48.7218 21.8452 22.21 12.14 12 1.0505 1.0313 1.01057 1.0076 1.021 0.992125 5 2.244 4.76 0.272 4.99 4.914 4.824 4.94 0.955 1.039 0.9 0.9537 0.9569];
%x = [44.303 18.564 10 10.1 12 1.1 1.0778 1.052 1.0574 1.0802 1.0803 4.3411 4.9527 4.2358 4.7605 4.0597 4.5901 4.1971 5 4.145 1 1.0125 1.025 1];
%x = [48.698 21.326 21.077 11.87 12 1.1 1.08 1.054 1.062 1.1 1.1 4.995 4.962 4.959 4.998 4.326 4.999 2.671 4.999 2.389 1.022 0.9 0.9665 0.9542];
%x = [51.97 15 10 10 12 1.0336 1.0113 0.9714 1.034 1.0993 1.1 4.984 5 5 5 4.65 5 5 5 5 1.1 1.053 1.07 1.065];
Qbus = [10 12 15 17 20 21 23 24 29];
Tbranch = [11 12 15 36];

data = loadcase(case30);
data.gen(2:6,2) = x(1:5);
data.gen(1:6,6) = x(6:11);
%data.bus(Qbus,4) = data.bus(Qbus,4)-x(12:20)';
data.bus(Qbus,6) = x(12:20);
data.branch(Tbranch,9) = x(21:24);

mpopt = mpoption('pf.enforce_q_lims',0,'verbose',0,'out.all',0);
result = runpf(data,mpopt);

rpowgen = [result.gen(1,2),x(1:5)];
%rpowgen = [198.55516,x(1:5)];
costcoeff = data.gencost(:,5:7);

fuelcost = sum(costcoeff(:,1)+costcoeff(:,2).*rpowgen'+costcoeff(:,3).*(rpowgen.^2)'); % be careful of sequence of coefficients

%Constraint finding 
Vmax = data.bus(:,12);
Vmin = data.bus(:,13);
genbus = data.gen(:,1);
pqbus = data.bus(:,1);
pqbus(genbus) = [];

Qmax = data.gen(:,4)/data.baseMVA;
Qmin = data.gen(:,5)/data.baseMVA;
QG = result.gen(:,3)/data.baseMVA;

PGSmax = data.gen(1,9);
PGSmin = data.gen(1,10);
PGS = result.gen(1,2);
PGSerr = (PGS<PGSmin)*(abs(PGSmin-PGS)/(PGSmax-PGSmin))+(PGS>PGSmax)*(abs(PGSmax-PGS)/(PGSmax-PGSmin));

blimit = data.branch(:,6);
Slimit = sqrt(result.branch(:,14).^2+result.branch(:,15).^2);
Serr = sum((Slimit>blimit).*abs(blimit-Slimit))/data.baseMVA;

% TO find the error in Qg of gen buses- inequality constraint
Qerr = sum((QG<Qmin).*(abs(Qmin-QG)./(Qmax-Qmin))+(QG>Qmax).*(abs(Qmax-QG)./(Qmax-Qmin)));
% TO find the error in V of load buses-inequality constraint
VI = result.bus(:,8);  %V of load buses-inequality constraint
VI_complx = VI.*(cosd(result.bus(:,9))+1i*sind(result.bus(:,9)));
vpvbus = VI_complx;
vpqbus = VI_complx;
vpvbus(pqbus) = [];
vpqbus(genbus) = [];
VI(genbus)=[];
Vmax(genbus)=[];
Vmin(genbus)=[];
VIerr = sum((VI<Vmin).*(abs(Vmin-VI)./(Vmax-Vmin))+(VI>Vmax).*(abs(Vmax-VI)./(Vmax-Vmin)));
VD = sum(abs(VI-1));

% Emission : gen_no. alpha beta gama omega miu d e
emcoeff = [
	1	0.04091 -0.05554 0.06490 0.000200 2.857 18 0.037;
	2	0.02543 -0.06047 0.05638 0.000500 3.333 16 0.038;
	3	0.04258 -0.05094 0.04586 0.000001 8.000 14 0.04;
    4	0.05326 -0.03550 0.03380 0.002000 2.000 12 0.045;
	5	0.04258 -0.05094 0.04586 0.000001 8.000 13 0.042;
	6	0.06131 -0.05555 0.05151 0.000010 6.667 13.5 0.041];

% VALVE EFFECT
%rpowgen12 = rpowgen(1:2);
%emcoeff12 = emcoeff(1:2,:);
%pgmin12 = data.gen(1:2,:);
%valveff = sum(abs(emcoeff12(:,7).*sin(emcoeff12(:,8).*(pgmin12(:,10)-rpowgen12')))); % if only 1st two have valve effects
valveff = sum(abs(emcoeff(:,7).*sin(emcoeff(:,8).*(data.gen(:,10)-rpowgen')))); % if all have valve effects

% BUS ADMITTANCE MATRICES
[Ybus,~,~] = makeYbus(data);
Ybuspq = Ybus;
Ybuspq(genbus,:) = [];
Ybuspvg = Ybuspq;
Ybuspq(:,genbus) = [];
Ybuspvg(:,pqbus) = [];
Fmat = -Ybuspq\Ybuspvg;

Lind = abs(1-(1./vpqbus).*(Fmat*vpvbus));
Lind_worst = max(Lind);

% OBJECTIVE FUNCTIONS
emission = sum(emcoeff(:,2)+emcoeff(:,3).*rpowgen'/100+emcoeff(:,4).*(rpowgen.^2/100^2)'...
    +emcoeff(:,5).*exp(emcoeff(:,6).*rpowgen'/100));
ploss = sum(result.branch(:,14)+result.branch(:,16));
qloss = sum(result.branch(:,15)+result.branch(:,17));

fuelvlvcost = fuelcost+valveff; % for CASE 3 only

error = [Qerr,VIerr,Serr,PGSerr];

objective_value = [PGS, fuelcost, fuelvlvcost, emission, ploss, VD, Lind_worst, qloss, VI'];

% f1 = fuelcost; % CASE 1: fuel cost only
%f1 = Lind_worst; % CASE 3: L-index % Case 2.4 in ref paper
%f1 = emission; % CASE 4: emission only 
%f1 = ploss; % CASE 5: power loss only
% f1 = fuelvlvcost; % CASE 6: gen with valve effect
% f1 = fuelcost+40*ploss; % CASE 7: fuel cost+power loss % Case 6 in ref
% f1 = fuelcost+100*VD; % CASE 8: fuel cost+voltage deviation, Case 7 in ref
% f1 = fuelcost+100*Lind_worst; % CASE 9: fuel cost+L-index, Case 8 in ref
f1 = fuelcost+19*emission+21*VD+22*ploss; % CASE 10: fuelcost+emission+VD+power loss
