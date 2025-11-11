 % The Ensemble of SF and SP constraint handling methods
clear all
global  nfeval %VD PGS emission ploss qloss Lind_worst VI fuelcost fuelvlvcost
format long e;
tic
D = 24;
NP = 50; 

Xmin = [20 15 10 10 12 0.95 0.95 0.95 0.95 0.95 0.95 0 0 0 0 0 0 0 0 0 0.9 0.9 0.9 0.9];
Xmax = [80 50 35 30 40 1.1 1.1 1.1 1.1 1.1 1.1 5 5 5 5 5 5 5 5 5 1.1 1.1 1.1 1.1];

Max_FES = 30000;
Max_Gen = 600; maxrun = 1; gn = 4; % gn is the no. of constraints
feval = []; contvar = zeros(maxrun,D);
pop = zeros(2*NP,D); val = zeros(); g = zeros(NP,gn); 
%pop2 = zeros(3*NP,D); val2 = zeros(); g2 = zeros(NP,gn);
Re = zeros(maxrun,3); W = zeros(Max_Gen,3); Splitfitness = zeros(Max_Gen,4);

Func = @pflow;

% JADE参数
mu_CR = 0.5;
mu_F = 0.5;
c = 0.1;
p = 0.05;
count_SA = 1;
SE = [];

 for runs=1:maxrun 
 
   fprintf('Run %d\n',runs);
   nfeval=0;
    
   pop=repmat(Xmin,NP,1)+repmat((Xmax-Xmin),NP,1).*rand(NP,D);
   val = zeros(NP,1); g = zeros(NP,gn); %PPB
   for i=1:NP
      [val(i,1),g(i,:)] = Func(pop(i,:));
   end
   
   nfeval=nfeval+NP;
 
 cons_max = [];

for gen=1:Max_Gen
    
    Sf = [];
    Scr = [];
       
    CR = [];
    F = [];
    
    for i=1:NP
        % 产生CR(i)
        CR(i)=normrnd(mu_CR,0.1);
        if CR(i)<0
            CR(i) = 0;
        end
        if CR(i)>1
            CR(i) = 1;
        end
        % 产生F(i)
        F(i)=cauchyrnd(mu_F,0.1);
        while F(i)<=0
            F(i)=cauchyrnd(mu_F,0.1);
        end
        
        if F(i)>1
            F(i) = 1;
        end
    end
    % CR排序
    [CR,~] = sort(CR);
    [~,indices] = sort(val);
    
     pid = ceil(NP*p);
     pbestid = randi(pid);
     Pbest = pop(indices(pbestid),:);
     for i =1:NP
        [Xr1,Xr2] = selectR(pop,i); 
        [m,~] = size(SE);
        if rand<=0.5 || m<1
            [Xr1,Xr2] = selectR(pop,i); 
            v(i,:) = pop(i,:)+F(i).*(Pbest-pop(i,:))+F(i).*(Xr1-Xr2);
        else
            indexSE = randi(m);
            v(i,:) = pop(i,:)+F(i).*(Pbest-pop(i,:))+F(i).*SE(indexSE,:);
        end
        % 交叉
        jrand = randi(D);
        crindex = find(indices==i);
        cr(i) = CR(crindex);
        for j=1:D
            if j==jrand || rand<cr(i)
                newpop(i,j) = v(i,j);
            else
                newpop(i,j) = pop(i,j);
            end
        end
     end
     
    newpop = boundConstraint(newpop,pop,Xmax,Xmin);
      
    newval = zeros(NP,1); newg = zeros(NP,gn); %PPB            
    for i=1:NP
      [newval(i,1),newg(i,:)] = Func(newpop(i,:));
    end             
    nfeval=nfeval+NP;
 
    A1 = vertcat(pop,newpop);%,newpop2);
    %A2 = vertcat(pop2,newpop2,newpop1);
    B1 = vertcat(val,newval);%,newval2);
    %B2 = vertcat(val2,newval2,newval1);
    G1 = vertcat(g,newg);%,newg2);
    %G2 = vertcat(g2,newg2,newg1);

   % Superiority of Feasible Solutions
    cons1=(G1>0).*G1; 
    cons_max=max([cons_max;cons1],[],1);
    nzindex1=find(cons_max~=0);
    
    if isempty(nzindex1)
        tcons1=zeros(size(A1,1),1);
    else
        tcons1=sum(cons1(:,nzindex1)./repmat(cons_max(:,nzindex1),size(A1,1),1),2)./sum(1./cons_max(:,nzindex1));
    end 
           
    ind11 = find(tcons1((NP+1):2*NP) < tcons1(1:NP));
    ind12 = find(tcons1((NP+1):2*NP) == 0);
    ind12 = ind12((tcons1(ind12) == 0));
    ind12 = ind12(B1(NP+ind12)<= B1(ind12));
    index1 = union(ind11,ind12);
    
    Len = length(index1);
    
    if Len>0
        for j=1:Len
            SE(count_SA,:) = A1(index1(j)+NP,:)-pop(index1(j),:);
            count_SA = count_SA+1;
            if count_SA>NP
                count_SA = 1;
            end
        end
     end
    
    pop(index1,:) = A1(NP+index1,:);
    val(index1,:) = B1(NP+index1,:);
    g(index1,:) = G1(NP+index1,:);
    
    Sf = F(index1);
    Scr = cr(index1);
    
    
       
    count = numel(Scr);
        
    if count<=1
        mu_CR = mu_CR;
        mu_F = mu_F;
    else
        mu_CR = (1-c)*mu_CR+c*mean(Scr);
        mu_F = (1-c)*mu_F+c*(sum(Sf.^2))/(sum(Sf));
    end

    %totalpop=[pop1;pop2];
    totalpop = pop;
    %totalval=[val1;val2];
    totalval = val;
    %totalg=[g1;g2];
    totalg = g;
     
    cons5=(totalg>0).*totalg; 
    cons_max=max([cons_max;cons5],[],1);
    nzindex5=find(cons_max~=0);
    
    if isempty(nzindex5)
        %tcons5=zeros(2*NP,1);
        tcons5=zeros(NP,1);
    else
        %tcons5=sum(cons5(:,nzindex5)./repmat(cons_max(:,nzindex5),2*NP,1),2)./sum(1./cons_max(:,nzindex5));
        tcons5=sum(cons5(:,nzindex5)./repmat(cons_max(:,nzindex5),NP,1),2)./sum(1./cons_max(:,nzindex5));
    end 
                
  feasindex=find(tcons5==0);
  if isempty(feasindex)
    [gbesttcons,ibest]=min(tcons5);
    gbestval=totalval(ibest);
    gbest = totalpop(ibest,:);
  else
    [gbestval,ibest]=min(totalval(feasindex));
    gbesttcons=tcons5(feasindex(ibest));
    gbest = totalpop(feasindex(ibest),:);
  end
  
  if(gen == 1)
      thegbestval = gbestval;
      thegbesttcons = gbesttcons;
  elseif((gbesttcons < thegbesttcons) || (gbesttcons==0 && thegbesttcons ==0 && gbestval < thegbestval))
      thegbestval = gbestval;
      thegbesttcons = gbesttcons;
  end

    W(gen,1)=thegbestval;
    W(gen,2)=thegbesttcons;
    W(gen,3)=nfeval;
%     Func(gbest); %only for multi-objectives
%     Splitfitness(gen,1) = nfeval; %only for multi-objectives
%     Splitfitness(gen,2) = thegbestval; %only for multi-objectives
%     Splitfitness(gen,3) = fuelcost; %only for multi-objectives
%     Splitfitness(gen,4) = Lind_worst; %only for multi-objectives
%     
    
 if(isempty(feval) && thegbesttcons == 0)
        feval = nfeval;
 end
%  thegbestval
%  fuelcost
%  Lind_worst
end

dlmwrite(strcat('W_',char(num2str(runs)),'.txt'),W,'precision','%0.4f','newline','pc');
Re(runs,1) = thegbestval;
Re(runs,2) = thegbesttcons;
Re(runs,3) = feval;

contvar(runs,:) = gbest;
 end
 
 pp = [min(Re(:,1)),max(Re(:,1)),mean(Re(:,1)),std(Re(:,1))];

dlmwrite('Result.txt',Re,'precision','%0.4f','newline','pc');
eval('save contvar');
toc

% % Print results
% [~,IND] = min(Re(:,1));
% disp('Control Variables');
% fprintf('\t %0.4f',contvar(IND,:));
% Func(contvar(IND,:));
% fprintf('\n Swing generator power %0.5f MW',PGS);
% fprintf('\n Cumulative voltage drop %0.5f p.u.',VD);
% fprintf('\n Emission %0.5f ton/h',emission);
% fprintf('\n Real power loss %0.4f MW',ploss);
% fprintf('\n Reative power loss %0.4f MVAr',qloss);
% fprintf('\n L-index %0.5f',Lind_worst);
% fprintf('\n Fuelcost %0.4f',fuelcost);
% fprintf('\n Fuelvlvcost %0.4f \n',fuelvlvcost);
% disp('Load bus voltage');
% fprintf('\t %0.5f',VI);