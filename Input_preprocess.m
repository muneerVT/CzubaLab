%% 
clear
clc
load C:\vt\research\Utah\Thesis\JD\dailyQ.mat; %synQ(:)=33;

%%
r_nodes=shaperead('C:\vt\research\Utah\Scott\v6\PLI_FinalV6\FullVB\Jordanelle\a009direct_input_nodes_snap_attV2.shp');
nodes = shaperead('C:\vt\research\Utah\Scott\v6\PLI_FinalV6\FullVB\Jordanelle\a009_debrisflow_pp_snap_netatt.shp');
%%
network = shaperead('C:\vt\research\Utah\Scott\v6\PLI_FinalV6\FullVB\Jordanelle\a009_network.shp');
LinkNum=length(network);
%%

GridID_orig = [network.GridID].';
GridID=(1:LinkNum)';

ToNode_orig = [network.ToLink].';

for i=1:LinkNum
if ToNode_orig(i) ~= 0
    ToNode(i)= GridID(find(ToNode_orig(i)==GridID_orig));
end 
end
ToNode=ToNode';

Length = [network.Length_m].';
mnelev=[network.dselev_m].';
mxelev=[network.uselev_m].';
Slope = [network.Slope].';
usarea_km=[network.usarea_km2].';
B=[network.riv_width].';
RES=[network.RES].';
OutletLinkID=GridID(ToNode_orig==0);

A_manual_change=find(GridID==ToNode);

%ToNode(OutletLinkID)=0; %manual
load('forSlp.mat')
% MANUAL CORRECTION (Seeing Network Shp) to ToNode
%% 
clear A_manual_change
Connect(1:LinkNum,1:LinkNum)=NaN;
Connect(1:LinkNum,1)=(1:LinkNum);
%%
for i=1:LinkNum
    j=1;
 while ToNode(Connect(i,j),1)~=ToNode(OutletLinkID)
           Connect(i,j+1)=find(ToNode(Connect(i,j),1)==GridID);
                j=j+1;
    end 
end
%Connect(Connect==0)=NaN;
%% IMPORT Q
%Date(1:6*365,1:3)=NaN; Qgage(1:6*365,1)=NaN; 

%% DF
NodeG_ID = [nodes.GridID].';

for i=1:length(NodeG_ID)
if NodeG_ID(i) ~= 0
    DF_LinkIn(i)= GridID(find(NodeG_ID(i)==GridID_orig));   
 end 
end
DF_LinkIn = DF_LinkIn';
DF_VolIn = [nodes.VolInput1].';

VolInDF=accumarray(DF_LinkIn,DF_VolIn);
%% Rusle
R_NodeG_ID = [r_nodes.GridID].';

for i=1:length(R_NodeG_ID)
if R_NodeG_ID(i) ~= 0
    RF_LinkIn(i)= GridID(find(R_NodeG_ID(i)==GridID_orig));    
 end 
end
RFs=40.6; % prev data
RF_LinkIn = RF_LinkIn';
RF_VolIn = 0.01*([r_nodes.R_HSD100].*RFs)';
RF_VolOut = sum(0.01*([r_nodes.R_HSD100].*(100-RFs))');
VolInRF=accumarray(RF_LinkIn,RF_VolIn);

%%
Agage= 230*2.58999 ; %km2 (Input)
%Q5(:)=12; % 7.34,10, 12, 13.72, 16.2, 17.8
Qgage=dailyQ;%Q_mixed; %if RP %if low: plus 1
%% if repeating
% Qgage(length(Qgage)+1:length(Qgage)*2,:)=Qgage; %to increase timesteps
% Qgage(length(Qgage)+1:length(Qgage)*2,:)=Qgage; %to increase timesteps
% Qgage(length(Qgage)+1:length(Qgage)*2,:)=Qgage;
%%
timesteps=length(Qgage);

a=0.84; %region 2, nearby site (input)

for t=1:timesteps
    for i=1:LinkNum
   Q(t,i)=Qgage(t)*(usarea_km(i)/Agage)^a;
         end
end
%% snyder (2013)

% % 2-yr flow
% 
% Qs=sort(Q,'descend');
% %Qs=Qs(~isnan(Qs));
% L=length(Qs);
% rank=(1:L)';
% pexceed=rank./(L+1).*100;
% 
% p=1/(365*2)*100;%2 yr
% Q2 = 10.^(interp1(pexceed,log10(Qs),p));

%% LD Q2

for i=1:LinkNum
   Q2(i)=66.4*(usarea_km(i)/Agage)^a; % from d/s gage Q2
end
%% %imp Q2 if req
dmean(1:LinkNum,1)=NaN;
for i=1:LinkNum
dmean(i,1)=0.045^(3/5)*Q2(1,i)^(3/5)*B(i,1)^(-3/5).*Slope(i,1)^(7/10)./0.05/1.65; %modified Q2(gage) n tauC* (input)
end
%modified Q2(gage) n tauC*
mean(dmean)
dmean(dmean<0.5*mean(dmean))=0.5*mean(dmean); mean(dmean)
clear  Qs L rank pexceed p a i t
%% stable high slope
% for i=1:LinkNum
%     if Slope(i)>0.02
%         dmean(i,1)=1.5*dmean(i,1);
%     end
% end
% 
% for i=1:LinkNum
%     if Slope(i)>0.04
%         dmean(i,1)=1.4*dmean(i,1);
%     end
% end

%% Distribution (Surface and Subsurface)

gsclass=8;
Dpsd=[0.5,2.83,5.66,16,45.3,90.5,181,362]'./1000;%m
%Fsfpsd=[0.046,0.024,0.048,0.274,0.278,0.166,0.1,0.064]'; % Twt surface
Fsfpsd=[0.0097,0.0041,0.015,0.0455,0.2053,0.3711,0.3034,0.0459]'; % RP surfac

% gravel calc
Gpsd=Dpsd(2:end);
for i=1:gsclass-1
Gsfpsd(i,1)= Fsfpsd(i+1)./(1-Fsfpsd(1));
end
Fsfpsdobs=Fsfpsd;
for i=1:length(Gpsd)
Gpsdphi(i)=log(Gpsd(i,1).*1000)/log(2); 
end
CumGsfpsd=cumsum(Gsfpsd);
% compute D16

for j=1:length(CumGsfpsd)
    if CumGsfpsd(j)<0.16 && CumGsfpsd(j+1)>=0.16
        D16=interp1(CumGsfpsd(j:j+1),Gpsdphi(j:j+1),0.16);
    end
end
for j=1:length(CumGsfpsd)
    if CumGsfpsd(j)<0.84 && CumGsfpsd(j+1)>=0.84
        D84=interp1(CumGsfpsd(j:j+1),Gpsdphi(j:j+1),0.84);
    end
end
stdphi=(D84-D16)/2; %making wider

for i=1:LinkNum
phimean(i,1)=log(dmean(i,1).*1000)/log(2); 
end


for i=1:LinkNum
pd(i,1)=makedist('normal','mu',(phimean(i,1)),'sigma',(stdphi));
end


x=Gpsdphi;
for i=1:LinkNum
    p(1:gsclass-1,i)=pdf(pd(i),x);%getting the values
end

P=sum(p);
for i=1:gsclass-1
    for j=1:LinkNum
        p(i,j)=p(i,j)./P(j); % total area=1
    end
end


s=0.046; %sand fraction (Input)
p=p';
Fsfpsd(LinkNum,gsclass)=NaN;
Fsfpsd(:,1)=s;

for i=1:LinkNum
    for j=2:gsclass
        Fsfpsd(i,j)=p(i,j-1).*(1-s);
    end
end

%%
Fsspsd=Fsfpsd; %subsurface

clear Gpsdphi CumGsfpsd Gsfpsd D84 D16 stdphi Gphimean GPhipsd Gpsd Gsfpsd j k i s P p pd phimean s

%%
%% Input
%timesteps=2192;

P_idx=cell(timesteps,LinkNum);
P_loc=cell(timesteps,LinkNum);
P_storage=cell(timesteps,LinkNum);
P_d=cell(timesteps,LinkNum);
P_vol=cell(timesteps,LinkNum);
P_tt=cell(timesteps,LinkNum);

theta= 0.25; %Assumed (dome above dmean)

VolIn=B.*Length.*theta;%m3


VolInd=repmat(VolIn.*4,1,gsclass).*Fsspsd;

pidx=1;
for i = 1:length(VolInd)
    for k = 1:gsclass
        %j=7-k+1;
        j=k;
        if VolInd(i,j) == 0
            continue
        end
        np=ceil(VolInd(i,j)./10);
        pvol=VolInd(i,j)./np;
        %P_loc{1,i}=cat(2,P_loc{1,i},repmat(0,1,np));
        P_loc{1,i}=cat(2,P_loc{1,i},rand(1,np));
        P_idx{1,i}=cat(2,P_idx{1,i},pidx+(0:1:np-1));
        pidx=pidx+np;
        P_storage{1,i}=cat(2,P_storage{1,i},zeros(1,np));%activate
        P_tt{1,i}=cat(2,P_tt{1,i},zeros(1,np));%all 0 sec
        P_vol{1,i}=cat(2,P_vol{1,i},repmat(pvol,1,np));
        P_d{1,i}=cat(2,P_d{1,i},repmat(Dpsd(j),1,np));
    end
end

pidxsf=pidx;

%surface
VolInd=repmat(VolIn.*3,1,gsclass).*Fsspsd;
for i = 1:length(VolInd)
    for k = 1:gsclass
        %j=7-k+1;
               j=k;
        if VolInd(i,j) == 0
            continue
        end
        np=ceil(VolInd(i,j)./10);
        pvol=VolInd(i,j)./np;
        %P_loc{1,i}=cat(2,P_loc{1,i},repmat(0,1,np));
        P_loc{1,i}=cat(2,P_loc{1,i},rand(1,np));
        P_idx{1,i}=cat(2,P_idx{1,i},pidx+(0:1:np-1));
        pidx=pidx+np;
        P_storage{1,i}=cat(2,P_storage{1,i},zeros(1,np));%activate
        P_tt{1,i}=cat(2,P_tt{1,i},zeros(1,np));%all 0 sec
        P_vol{1,i}=cat(2,P_vol{1,i},repmat(pvol,1,np));
        P_d{1,i}=cat(2,P_d{1,i},repmat(Dpsd(j),1,np));
    end
end
%% sedpile/US supply
% for i=1:LinkNum
%  usid=find(Connect(:,2)==i);
% end

%select GSD

% Fdfpsd=[0.25,0.06225,0.4168,0.187866667,0.048616667,0.034466667,0,0]'; % debris flow A - Twitchell
%Fdfpsd=[0.78,0.019366667,0.1317,0.068933333,0,0,0,0]'; % debris flow B - Yellowstone_1988
Fdfpsd=[0.14,0.0332,0.228933333,0.161383333,0.168616667,0.15585,0.112016667,0]'; % debris flow C - Coal_Hollow
%Fdfpsd=[0.53,0.02075,0.140316667,0.082766667,0.094466667,0.062766667,0.068933333,0]'; % debris flow D -Trailmountain
Fdfpsd=repmat(Fdfpsd',LinkNum,1);

%Fdfpsd=Fsfpsd;
%Fdfpsd=repmat(Fdfpsd,LinkNum,1);

pidxus=pidx;

maxparcelsize=10;
% QQss=sum(Qss,1);
% VolInd=QQss.*3;

VolIn=B.*Length.*theta;%m3
% 
% 
USlink = setdiff(GridID,ToNode);
% %VolInd=repmat(VolIn.*2*timesteps/length(USlink),1,gsclass).*Fdfpsd;
VolInd=repmat(VolIn.*0.3.*B*timesteps/length(USlink),1,gsclass).*Fdfpsd; %scale to cap

% QQss=sum(Qss,1);
% VolInd=8000.*QQss;
% VolInd = repmat(VolInd,LinkNum,1);

% VolInd=VolIn(202)*40.*Fdfpsd(202,:); %u/s
% VolInd(202,1:8)=VolInd;
for i=1:LinkNum
    if ~ismember(i,USlink,'legacy')
        VolInd(i,1:gsclass)=0;
    end
end


for i = 1:LinkNum %:length(VolInd)
    for k = 1:gsclass
        %j=7-k+1;
        j=k;
        if VolInd(i,j) == 0
            continue
        end
        np=ceil(VolInd(i,j)./maxparcelsize);
        pvol=VolInd(i,j)./np;
        %P_loc{1,i}=cat(2,P_loc{1,i},repmat(0,1,np));
        P_loc{1,i}=cat(2,P_loc{1,i},rand(1,np));
        P_idx{1,i}=cat(2,P_idx{1,i},pidx+(0:1:np-1));
        pidx=pidx+np;
        P_storage{1,i}=cat(2,P_storage{1,i},zeros(1,np));%activate
        P_tt{1,i}=cat(2,P_tt{1,i},zeros(1,np));%all 0 sec
        P_vol{1,i}=cat(2,P_vol{1,i},repmat(pvol,1,np));
        P_d{1,i}=cat(2,P_d{1,i},repmat(Dpsd(j),1,np));
    end
end
%% RUSLE input

Frfpsd=[1,0,0,0,0,0,0,0]'; % sand fraction only
Frfpsd=repmat(Frfpsd',LinkNum,1);

pidxrf=pidx;

VolInd=repmat(VolInRF,1,gsclass).*Frfpsd; %prepare VolInRF
for i = 1:length(VolInd)
    for k = 1:gsclass
        %j=7-k+1;
        j=k;
        if VolInd(i,j) == 0
            continue
        end
        np=ceil(VolInd(i,j)./10);
        pvol=VolInd(i,j)./np;
        %P_loc{1,i}=cat(2,P_loc{1,i},repmat(0,1,np));
        P_loc{1,i}=cat(2,P_loc{1,i},rand(1,np));
        P_idx{1,i}=cat(2,P_idx{1,i},pidx+(0:1:np-1));
        pidx=pidx+np;
        P_storage{1,i}=cat(2,P_storage{1,i},zeros(1,np));%activate
        P_tt{1,i}=cat(2,P_tt{1,i},zeros(1,np));%all 0 sec
        P_vol{1,i}=cat(2,P_vol{1,i},repmat(pvol,1,np));
        P_d{1,i}=cat(2,P_d{1,i},repmat(Dpsd(j),1,np));
    end
end
%% Debris flow

%select GSD
Fdfpsd=[0.25,0.06225,0.4168,0.187866667,0.048616667,0.034466667,0,0]'; % debris flow A - Twitchell
%Fdfpsd=[0.78,0.019366667,0.1317,0.068933333,0,0,0,0]'; % debris flow B - Yellowstone_1988
%Fdfpsd=[0.14,0.0332,0.228933333,0.161383333,0.168616667,0.15585,0.112016667,0]'; % debris flow C - Coal_Hollow
%Fdfpsd=[0.53,0.02075,0.140316667,0.082766667,0.094466667,0.062766667,0.068933333,0]'; % debris flow D -Trailmountain

Fdfpsd=repmat(Fdfpsd',LinkNum,1);

pidxdf=pidx;

VolInd=repmat(VolInDF,1,gsclass).*Fdfpsd;

for i = 1:length(VolInd)
    for k = 1:gsclass
        %j=7-k+1;
        j=k;
        if VolInd(i,j) == 0
            continue
        end
        np=ceil(VolInd(i,j)./10);
        pvol=VolInd(i,j)./np;
        %P_loc{1,i}=cat(2,P_loc{1,i},repmat(0,1,np));
        P_loc{1,i}=cat(2,P_loc{1,i},rand(1,np));
        P_idx{1,i}=cat(2,P_idx{1,i},pidx+(0:1:np-1));
        pidx=pidx+np;
        P_storage{1,i}=cat(2,P_storage{1,i},zeros(1,np));%activate
        P_tt{1,i}=cat(2,P_tt{1,i},zeros(1,np));%all 0 sec
        P_vol{1,i}=cat(2,P_vol{1,i},repmat(pvol,1,np));
        P_d{1,i}=cat(2,P_d{1,i},repmat(Dpsd(j),1,np));
    end
end

%% saving Inputs

save JD_dailyQ_dfA.mat
%%
%Q=2*Q;
%save RB_2Q_C_100G_inputV6.mat

