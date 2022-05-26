%% Initialize variables

%%
clear all
close all
clc
%% Load
load JD_dailyQ_dfA_ths.mat
theta=0.25; %Q=Q.*1.4;
%%
daystp=1;%182.5;%18.25; %number of days per timestep
dt = 1*24*60*60*daystp; %seconds in daystp number of days

tmax=dt*(timesteps-1)/60/60/24/365; %years, max time of simulation
time=(0:daystp/365:tmax)'; %years, time of current step

%theta=0.25;%0.064;%0.1;%2D90, 64mm, units in meters
g=9.81;%m/s2
R=1.65;
rho=1000;%kg/m3

OutVol(1:timesteps,1)=0;

bH=(.4655*usarea_km.^0.3549)/3.28;

%compute capacity for all links
capacity=B.*Length.*theta;%m3

for i=1:LinkNum
    index=Connect(i,~isnan(Connect(i,:)));
    for j=2:1:length(index)
       
        if capacity(index(j))<(0.7*capacity(index(j-1)))  %setting main condition
        capacity(index(j))=0.7*capacity(index(j-1));
        else
            capacity(index(j))= capacity(index(j)); 
        end
        capacity(index(1))=capacity(index(1));
    end
    % what about at confluence
end

USlink = setdiff(GridID,ToNode);
while capacity(USlink)<mean(capacity)
    capacity(USlink)=mean(capacity);
end


Slope_raw=Slope;
Slope_raw(Slope_raw<0.001)=0.001;
%modify elevations so they are consistent with adjusted slopes
% mxelevmod(1:LinkNum,1)=NaN;
% for i = 1:LinkNum
%     idx=fliplr(Connect(i,~isnan(Connect(i,:))));
%     for j=1:length(idx)
%         if idx(j)==OutletLinkID
%             mxelevmod(idx(j),1)=mxelev(OutletLinkID)+Slope_cond(idx(j)).*Length(idx(j));
%         else
%             mxelevmod(idx(j),1)=mxelevmod(idx(j-1),1)+Slope_cond(idx(j)).*Length(idx(j));
%         end
%     end
%     clear idx
% end
% clear i j

Elev=repmat(mxelev',timesteps,1);
% %Slope - slope at current timestep
% %slope - original slope
% slope=Slope;
% %Slope=repmat(slope',timesteps,1);

Dg_Lp(1:timesteps,1:LinkNum)=NaN;
Lp(1:timesteps,1:LinkNum)=NaN;
%Lp=0.25;%0.21;%0.4; %porosity
%Lp(1,:)=[0.13+0.21./(((dmean.*1000)+0.002).^0.21)]';%Wu and Wang 2006

Lake(1:LinkNum,1)=0;

abrasion_rate = 0; %0.02 ./ 1000; %km-1 to m-1 
lnkvol(1:timesteps,1:LinkNum)=NaN;
Dg(1:timesteps,1:LinkNum)=NaN;

Dcsd=[2,4,8,32,64,128,256,512]'./1000;%m
DcsdG=Dcsd(2:end);


H(1:timesteps,1:LinkNum)=NaN;
Dg(1:timesteps,1:LinkNum)=NaN;
DG(1:timesteps,1:LinkNum)=NaN;
D84(1:timesteps,1:LinkNum)=NaN;
D50(1:timesteps,1:LinkNum)=NaN;
D16(1:timesteps,1:LinkNum)=NaN;
D84G(1:timesteps,1:LinkNum)=NaN;
D50G(1:timesteps,1:LinkNum)=NaN;
D16G(1:timesteps,1:LinkNum)=NaN;
% Rickemann and Recking SP
a1=6.5;
a2=2.5;
So=Slope;

%% random sort
% 
for i=1:LinkNum
% sub
subidx=find(P_idx{1,i}<pidxsf);
n = length(subidx);
if n ~= 0
    ridx=randperm(n); 
    P_idx{1,i}(subidx)=P_idx{1,i}(subidx(ridx));
    P_vol{1,i}(subidx)=P_vol{1,i}(subidx(ridx));
    P_d{1,i}(subidx)=P_d{1,i}(subidx(ridx)); 
end

%sf
sfidx=find(and(P_idx{1,i}<pidxus,P_idx{1,i}>=pidxsf));
n = length(sfidx);
if n ~= 0
    ridx=randperm(n); 
    P_idx{1,i}(sfidx)=P_idx{1,i}(sfidx(ridx));
    P_vol{1,i}(sfidx)=P_vol{1,i}(sfidx(ridx));
    P_d{1,i}(sfidx)=P_d{1,i}(sfidx(ridx)); 
end

%% us
usidx=find(and(P_idx{1,i}<pidxdf,P_idx{1,i}>=pidxus));
n = length(usidx);
if n ~= 0
    ridx=randperm(n); 
    P_idx{1,i}(usidx)=P_idx{1,i}(usidx(ridx));
    P_vol{1,i}(usidx)=P_vol{1,i}(usidx(ridx));
    P_d{1,i}(usidx)=P_d{1,i}(usidx(ridx)); 
end
%% df


dfidx=find(P_idx{1,i}>=pidxdf);
n = length(dfidx);
if n ~= 0
    ridx=randperm(n); 
    P_idx{1,i}(dfidx)=P_idx{1,i}(dfidx(ridx));
    P_vol{1,i}(dfidx)=P_vol{1,i}(dfidx(ridx));
    P_d{1,i}(dfidx)=P_d{1,i}(dfidx(ridx)); 
end
end

%% Map input structure
% so inputs can always be added on top of what is in the system

P_idx_IN=P_idx;
P_loc_IN=P_loc;
P_storage_IN=P_storage;
P_vol_IN=P_vol;
P_d_IN=P_d;
P_tt_IN=P_tt;

P_idx=cell(timesteps,LinkNum);
P_loc=cell(timesteps,LinkNum);
P_storage=cell(timesteps,LinkNum);
P_vol=cell(timesteps,LinkNum);
P_d=cell(timesteps,LinkNum);
P_tt=cell(timesteps,LinkNum);

pidxcap=pidx;

%% Run Model

tic

for t = 1:timesteps-1 %step through each timestep
    t 
    % Generate autogenic inputs here? If desired.
    
%     % Add new inputs at given timestep on top of what is already in system
%     % background inputs
    for i=1:LinkNum
        P_loc{t,i}=cat(2,P_loc{t,i},P_loc_IN{t,i});
        P_idx{t,i}=cat(2,P_idx{t,i},P_idx_IN{t,i});
        P_storage{t,i}=cat(2,P_storage{t,i},P_storage_IN{t,i});
        P_vol{t,i}=cat(2,P_vol{t,i},P_vol_IN{t,i});
        P_d{t,i}=cat(2,P_d{t,i},P_d_IN{t,i});
        P_tt{t,i}=cat(2,P_tt{t,i},P_tt_IN{t,i});      
    end
    
%     for i=1
%        P_vol{t+1,i}=cat(2,P_vol{t,i},P_vol_IN{t,i});
%     end
        
    
    
%     %input pulse here
%     if t==2001 %add pulse at 100 yrs
%         for i=1:LinkNum
%             P_loc{t,i}=cat(2,P_loc{t,i},P_loc_PIN{t,i});
%             P_idx{t,i}=cat(2,P_idx{t,i},P_idx_PIN{t,i});
%             P_storage{t,i}=cat(2,P_storage{t,i},P_storage_PIN{t,i});
%             P_vol{t,i}=cat(2,P_vol{t,i},P_vol_PIN{t,i});
%         end
%         clear P_loc_PIN P_idx_PIN P_storage_PIN P_vol_PIN
%     end
    
    % Determine capacity of links, active parcels, and adjust bed slopes
    % based on storage volume
    CapacitySlope_Utah_2222022S

    %Velocity=0.05./sqrt(9.81)./1.65./1.65./D.*U.^2.*H.^(1/2).*Slope.^(3/2)./theta.*0.175;%m/s, realtime
    %STTime=Length./Velocity;%seconds, travel time through each link
    %STTime -> P_tt{t,i}(p)

    for i = 1:LinkNum %step through each link
        
        % if current link is empty go to the next link
        if isempty(P_loc{t,i})
            continue
        end
        
        if RES(i)==1
            continue
        end
        
%         if Lake(i) % if a lake, do not move sediment
%             continue
%         end
        
        %capacity and slope here? no, needs to be above i loop because
        %slopes are updated from more than just its own link
        
        for p = 1:length(P_loc{t,i})
            
            % if parcel is inactive then update and go to the next parcel
            if P_storage{t,i}(1,p)
                P_loc{t+1,i}=cat(2,P_loc{t+1,i},P_loc{t,i}(1,p));%keep loc
%                 if P_loc{t,i}(1,p)<0
%                     1;
%                 end
                P_idx{t+1,i}=cat(2,P_idx{t+1,i},P_idx{t,i}(1,p));%keep idx
                P_storage{t+1,i}=cat(2,P_storage{t+1,i},0);% make active as 
                %this is determined for each timestep
                
                P_vol{t+1,i}=cat(2,P_vol{t+1,i},P_vol{t,i}(1,p));%keep vol
                P_d{t+1,i}=cat(2,P_d{t+1,i},P_d{t,i}(1,p));%keep d
                P_tt{t+1,i}=cat(2,P_tt{t+1,i},0);% make 0 as 
                %this is determined for each timestep
                continue
            end
            
            % Enter long term storage here?
            
            % Move parcels downstream
            sm=[];%time to move through current and downstream links
            csm=[];%cumulative time to move through current and ds links
            sm=P_tt{t,i}(p)*(1-P_loc{t,i}(1,p));%seconds, time to move out of current link
            csm=sm;%seconds
            ii=i;
            while csm(end)<=dt %loop until cumulative time toward outlet exceeds timestep
                ii=Connect(ii,2); %set ds link
                if isnan(ii) %if downstream is the outlet
                    %leave outlet before dt ends
                    sm=cat(1,sm,NaN);
                    csm=cat(1,csm,sum(sm));
                    OutVol(t,1)=OutVol(t,1)+P_vol{t,i}(1,p);
                    break
                end
                
                if RES(ii)==1 %if downstream is the reservoir
                    %leave system before dt ends
                    sm=cat(1,sm,NaN);
                    csm=cat(1,csm,sum(sm));
                    OutVol(t,1)=OutVol(t,1)+P_vol{t,i}(1,p);
                    break
                end
                
%                 if Lake(ii) %if move into a lake
%                     sm=cat(1,sm,NaN);
%                     csm=cat(1,csm,sum(sm));                   
%                    break 
%                 end
                %sm=cat(1,sm,STTime(ii,1));%add time through next link
                %EDIT IN THE FUTURE
                sm=cat(1,sm,P_tt{t,i}(p)./Length(i).*Length(ii));%travel at same velocity but in ds link
                csm=cat(1,csm,sum(sm));%cumulative time to move through all subsequent links
            end
            % update parcel location
            if ~isnan(csm(end)) %check to make sure parcel is still in the system
                pi=Connect(i,length(csm)); %parcel link
                if length(csm)==1 %still in same link
                    pl=(P_tt{t,i}(p)*P_loc{t,i}(1,p)+dt)/P_tt{t,i}(p);%update location from current ploc
                else %moved to a ds link
                    pl=1-((csm(end)-dt)/P_tt{t,i}(p));%update location from beginning of link
                    if pl<1 %overshooting the next link, UPDATE IN FUTURE
                        pl=1;
                    end
                    %note csm(end)-dt computes time remaining to move
                    %through the rest of the link
                end
                              
                P_loc{t+1,pi}=cat(2,P_loc{t+1,pi},pl);
%                 if pl<0
%                     1;
%                 end
                P_idx{t+1,pi}=cat(2,P_idx{t+1,pi},P_idx{t,i}(1,p));
                P_storage{t+1,pi}=cat(2,P_storage{t+1,pi},0);%activate
                
                distance_traveled = P_tt{t,i}(p).*dt; %m
                new_vol = P_vol{t,i}(1,p).*exp(distance_traveled.*(-abrasion_rate)); %m3
                new_d = (P_d{t,i}(1,p)) *(new_vol/P_vol{t,i}(1,p)).^(1/3); %m
                
                if P_d{t+1,pi}>0.002
                    P_vol{t+1,pi}=cat(2,P_vol{t+1,pi},new_vol);
                    P_d{t+1,pi}=cat(2,P_d{t+1,pi},new_d);
                else
                    P_vol{t+1,pi}=cat(2,P_vol{t+1,pi},P_vol{t,i}(1,p));
                    P_d{t+1,pi}=cat(2,P_d{t+1,pi},P_d{t,i}(1,p));
                end


                P_tt{t+1,pi}=cat(2,P_tt{t+1,pi},0);%set to 0
                
                
%                                 if i==1 && pi>1
%                     % Keep US at capacity
%                     % any parcels that leave US link, go back into US link
%                     P_loc{t+1,i}=cat(2,P_loc{t+1,i},0);
%                     P_idx{t+1,i}=cat(2,P_idx{t+1,i},pidx);
%                     P_storage{t+1,i}=cat(2,P_storage{t+1,i},0);%activate
%                     P_vol{t+1,i}=cat(2,P_vol{t+1,i},P_vol{t,i}(1,p));
%                     P_d{t+1,i}=cat(2,P_d{t+1,i},P_d{t,i}(1,p));
%                     P_tt{t+1,i}=cat(2,P_tt{t+1,i},0);%set to 0
%                     pidx=pidx+1;
%                 end
                
                
%                 if i==1 && pi>1
%                     % Keep US at capacity
%                     % any parcels that leave US link, go back into US link
%                     P_loc{t+1,i}=cat(2,P_loc{t+1,i},0);
%                     P_idx{t+1,i}=cat(2,P_idx{t+1,i},pidx);
%                     P_storage{t+1,i}=cat(2,P_storage{t+1,i},0);%activate
% 
%                     distance_traveled = P_tt{t,i}(p).*dt; %m
%                     new_vol = P_vol{t,i}(1,p).*exp(distance_traveled.*(-abrasion_rate)); %m3
%                     new_d = (P_d{t,i}(1,p)) *(new_vol/P_vol{t,i}(1,p)).^(1/3); %m
%                     P_vol{t+1,pi}=cat(2,P_vol{t+1,pi},new_vol);
%                     P_d{t+1,pi}=cat(2,P_d{t+1,pi},new_d);
% %                     P_vol{t+1,i}=cat(2,P_vol{t+1,i},P_vol{t,i}(1,p));
% %                     P_d{t+1,i}=cat(2,P_d{t+1,i},P_d{t,i}(1,p));
%                     P_tt{t+1,i}=cat(2,P_tt{t+1,i},0);%set to 0
%                     pidx=pidx+1;
%                 end
%                 
%                P_mass{t+1,pi}=cat(2,P_mass{t+1,pi},P_mass{t,i}(1,p).*exp(-kd*dt));
            end
%             %assign arrival times
%             for ai=1:length(csm)-1 %as long the parcel moved to a new link
%                 if isnan(Connect(i,ai+1)) %left the system
%                     OutArrival=cat(1,OutArrival,...
%                         time(t,1)+csm(ai,1)/60/60/24/365);
%                 else %arrival times at each link
%                     L_arrival{Connect(i,ai+1),1}=cat(2,L_arrival{Connect(i,ai+1),1},...
%                         time(t,1)+csm(ai,1)/60/60/24/365);
%                 end
%             end
            
        end %p loop
    end %i loop 
    
%     if t==2001
%         lnkvol=cellfun(@sum,P_vol); %sum the volume in each link
%         numpar=cellfun(@length,P_loc); %compute number of parcel in each link
%         seddepth=lnkvol./repmat(Length',timesteps,1)./repmat(B',timesteps,1);%m
%         save('BE_NHD_MartinLakes4_BRU_baseline_100yr_v1.mat','seddepth','lnkvol','Elev','OutVol','OutArrival','L_arrival','-v7.3');
%         clear seddepth
%     end

% handle overshooting, placing P_loc=1 to P_loc=0 in next link
for i=1:LinkNum           
            if sum(P_loc{t+1,i}==1)
                
                P_idx{t+1,Connect(i,2)}=cat(2,P_idx{t+1,Connect(i,2)},P_idx{t+1,i}(P_loc{t+1,i}==1));
                P_storage{t+1,Connect(i,2)}=cat(2,P_storage{t+1,Connect(i,2)},P_storage{t+1,i}(P_loc{t+1,i}==1));
                P_vol{t+1,Connect(i,2)}=cat(2,P_vol{t+1,Connect(i,2)},P_vol{t+1,i}(P_loc{t+1,i}==1));
                P_d{t+1,Connect(i,2)}=cat(2,P_d{t+1,Connect(i,2)},P_d{t+1,i}(P_loc{t+1,i}==1));
                P_tt{t+1,Connect(i,2)}=cat(2,P_tt{t+1,Connect(i,2)},P_tt{t+1,i}(P_loc{t+1,i}==1));
                P_loc{t+1,Connect(i,2)}=cat(2,P_loc{t+1,Connect(i,2)},P_loc{t+1,i}(P_loc{t+1,i}==1).*0);
                
                P_idx{t+1,i}(P_loc{t+1,i}==1)=[];
                P_storage{t+1,i}(P_loc{t+1,i}==1)=[];
                P_vol{t+1,i}(P_loc{t+1,i}==1)=[];
                P_d{t+1,i}(P_loc{t+1,i}==1)=[];
                P_tt{t+1,i}(P_loc{t+1,i}==1)=[];
                P_loc{t+1,i}(P_loc{t+1,i}==1)=[];
            end
end
lnkvol(t,:)=cellfun(@sum,P_vol(t,:));
% 
% for i=1:LinkNum
%     Ppulse(t,i)=sum(P_idx{t,i}>=ppidxstart);
%     Ppvol(t,i)=sum(P_vol{t,i}(P_idx{t,i}>=ppidxstart));
% end
% 
% %clear contents to conserve space
% for i=1:LinkNum
%     P_idx{t,i}=[];
%     P_loc{t,i}=[];
%     P_storage{t,i}=[];
%     P_vol{t,i}=[];
%     
%     P_idx_IN{t,i}=[];
%     P_loc_IN{t,i}=[];
%     P_storage_IN{t,i}=[];
%     P_vol_IN{t,i}=[];
% end

end %t loop


clear sm csm ii pi pl ai

toc

%OutVol(1)=OutVol(1)+ RF_VolOut + sum(VolInRF (RES==1))+ sum(VolInDF (RES==1));
ResVol=cumsum(OutVol);
%%
%
% lnkmass=cellfun(@sum,P_mass); %sum the mass in each link
% lnkconc=lnkmass./(repmat((Length.*B.*H)',timesteps,1)).*1000;%mg/L
% numpar=cellfun(@length,P_loc); %compute number of parcel in each link

%
%lnkvol=cellfun(@sum,P_vol); %sum the volume in each link
%numpar=cellfun(@length,P_loc); %compute number of parcel in each link
%B=(0.0238).*(usarea).^(0.3397);
%numparconc=numpar./repmat(Length',timesteps,1)./repmat(B',timesteps,1).*1000;
t=timesteps;
lnkvol(t,:)=cellfun(@sum,P_vol(t,:));
for i=1:LinkNum
    if ~isempty(P_vol{t,i})
        Dg_Lp(t,i)=(2.^(sum(P_vol{t,i}./sum(P_vol{t,i}).*...
            (log10(P_d{t,i})./log10(2)))));
        Lp(t,i)=0.13+0.21./(((Dg_Lp(t,i).*1000)+0.002).^0.21);%Wu and Wang 2006
    end
end
seddepth=lnkvol./repmat(Length',timesteps,1)./repmat(B',timesteps,1)./(1-Lp);%m

%save MethowOutput.mat

%sum(sum(cellfun(@sum,P_storage),2))
%
%save('TusharOutput.mat','lnkvol','Dg','P_loc','P_vol','P_d','-v7.3');

%save('TusharOutput14.mat','-v7.3');

%% Compute sediment depth from DF parcels
lnkvoldf(1:timesteps,1:LinkNum)=NaN;
for i=1:LinkNum
    for tt=1:timesteps
        lnkvoldf(tt,i)=sum(P_vol{tt,i}(P_idx{tt,i}>=pidxdf));
    end
end
seddepthdf=lnkvoldf./repmat(Length',timesteps,1)./repmat(B',timesteps,1)./(1-Lp);%m

%% Compute surface sand fraction
Fs(1:timesteps,1:LinkNum)=NaN;
for t=1:timesteps
    for i=1:LinkNum
        if isempty(P_loc{t,i})
            continue
        end
        actidx=P_storage{t,i}==0;
        actvol=sum(P_vol{t,i}(actidx));
        
        actsandidx=and(P_storage{t,i}==0,P_d{t,i}==Dpsd(1));
        actsandvol=sum(P_vol{t,i}(actsandidx));
        Fs(t,i)=actsandvol./actvol;
    end
end

%%
%% TIMESERIES PLOTS
%%
%%
plot((1:timesteps)./365,ResVol,'k','LineWidth',2)
xlim([0 timesteps./365])
%ylim([0,12000])
ylabel('Total Volume within reservoir (no porosity), m^3')
xlabel('Time since debris flow input, years')

lnkvolgs(1:timesteps,1:LinkNum,1:gsclass)=NaN;
for j=1:gsclass
    for i=1:LinkNum
        for tt=1:timesteps
            %only works when no abrasion! otherwise need to bound P_d
            lnkvolgs(tt,i,j)=sum(P_vol{tt,i}(P_d{tt,i}==Dpsd(j)));
        end
    end
end

%mineral volume, not including porosity
volinbasings=sum(lnkvolgs,2);
totalvolinbasings=sum(volinbasings,3);

figure; hold on; box on
for i=1:gsclass-1
plot((1:timesteps)./365,volinbasings(:,1,i))
end
ylabel('Volume of bed input within network (no porosity), m^3')
xlabel('Time since flow input, years')
plot((1:timesteps)./365,totalvolinbasings(:,1,1),'k','LineWidth',2)
legend(cat(1,num2str(Dpsd(1:end-1).*1000),' all'))
%ylim([0,1])
xlim([0 timesteps./365])
%set(gca,'YScale','log')
%% Plot seddepth vs. t for all links
figure; hold on; box on
plot((1:timesteps)./365,seddepth)
ylabel('Total sediment depth, m^3')
xlabel('Time since flow input, years')
xlim([0 timesteps./365])

%% Plot seddepthdf vs. t for all links
figure; hold on; box on
plot((1:timesteps)./365,seddepthdf)
ylabel('Accumulation from debris flow sediment, m^3')
xlabel('Time since flow input, years')
xlim([0 timesteps./365])

%% Plot Dg vs. t for all links
% includes sand fraction in calculation
figure; hold on; box on
plot((1:timesteps)./365,Dg)
set(gca,'YScale','log')
ylabel('Dg (includes sand), m')
xlabel('Time since flow input, years')
xlim([0 timesteps./365])

%% Plot DG vs. t for all links
% excludes sand fraction in calculation
figure; hold on; box on
plot((1:timesteps)./365,DG)
set(gca,'YScale','log')
ylabel('DG (excludes sand), m')
xlabel('Time since flow input, years')
xlim([0 timesteps./365])

%% Plot Fs vs. t for all links
% surface sand fraction
figure; hold on; box on
plot((1:timesteps)./365,Fs)
ylabel('Surface sand fraction (Fs), m')
xlabel('Time since flow input, years')
xlim([0 timesteps./365])

%%
%% LINKSERIES
%%
%% Plot seddepth for all links
figure; hold on; box on
plot(seddepth')
ylabel('Total sediment depth, m^3')
xlabel('Link index')
xlim([1 LinkNum])

%% Plot seddepthdf for all links
figure; hold on; box on
plot(seddepthdf')
ylabel('Accumulation from debris flow sediment, m^3')
xlabel('Link index')
xlim([1 LinkNum])

%% Plot Dg for all links
% includes sand fraction in calculation
figure; hold on; box on
plot(Dg')
set(gca,'YScale','log')
ylabel('Dg (includes sand), m')
xlabel('Link index')
xlim([1 LinkNum])

%% Plot DG for all links
% excludes sand fraction in calculation
figure; hold on; box on
plot(DG')
set(gca,'YScale','log')
ylabel('DG (excludes sand), m')
xlabel('Link index')
xlim([1 LinkNum])

%% Plot Fs for all links
% surface sand fraction
figure; hold on; box on
plot(Fs')
ylabel('Surface sand fraction (Fs), m')
xlabel('Link index')
xlim([1 LinkNum])


% %%
% %% Network
% %%
% %% Determine incremental length spatially though links 2-D
% %length through link
% for i=1:LinkNum
%     network(i).L=network(i).X;
% end
% for i=1:LinkNum
%     network(i).L(1)=0;
%     for j=2:length(network(i).X)-1
%         network(i).L(j)=sqrt((network(i).X(j)-network(i).X(j-1)).^2+...
%             (network(i).Y(j)-network(i).Y(j-1)).^2);
%     end
% end
% for i=1:LinkNum
%     for j=3:length(network(i).X)-1
%         network(i).L(j)=network(i).L(j)+network(i).L(j-1);
%     end
%     network(i).L=network(i).L./nanmax(network(i).L);
% end
% 
% %% Load additional shapefiles
% boundary = shaperead('C:\Users\jczuba\Documents\Projects\UtahWildfire\GIS\analysis\Watersheds\080\a080_watershed2.shp');
% reservoir = shaperead('C:\Users\jczuba\Documents\Projects\UtahWildfire\GIS\analysis\Watersheds\080\a080_reservoir.shp');
% 
% %% Plot parcels on the network
% % configure file for each network by changing t, anchor, xlegscale, and
% % ylegscale
% % anchor is the relative position of the upper left corner of the
% % bounding box of the entire legend to the relative position of the
% % bounding box of the watershed boundary (starting in the lower left
% % corner). [1 1] would put the legend 1 watershed bounding box up and to
% % the right.
% % xlegscale and ylegscale are scale factors to help spread out the legend
% % entries if they get too crowded.
% t=80;%1+365*0;
% %anchor=[0.6,0.3];xlegscale=1.0;ylegscale=1;%001
% %anchor=[0.75,0.5];xlegscale=1.4;ylegscale=1;%076
% %anchor=[-0.2,0.3];xlegscale=1;ylegscale=1;%074
% anchor=[0.5,0.5];xlegscale=2.2;ylegscale=1;%080
% %anchor=[-0.1,0.15];xlegscale=1;ylegscale=1.4;%071
% %anchor=[-0.4,0.3];xlegscale=1;ylegscale=1;%061
% 
% BMS_PlotParcelsUtah
% 
% %% Create a movie of parcels moving on the network
% vfilename='a080_initialrun_6yr';
% % anchor, xlegscale, and ylegscale should be created and set in block above
% BMS_MovieParcelsUtah
% 
% %%
% %% PATHWAY
% %%
% %% Compute distance from every link to the outlet
% Dist(1:LinkNum,1)=NaN;
% for i=1:LinkNum
%     lastConn=find(Connect(i,:)==OutletLinkID);
%     Dist(i,1)=sum(Length(Connect(i,1:lastConn)),1);
% end
% 
% %% Elevation profile -- all links
% figure; hold on; box on
% for i=1:LinkNum
%     index=Connect(i,~isnan(Connect(i,:)));
%     plot([Dist(index,1)./1000; 0],...
%         [mxelev(index,1); mnelev(index(end),1)],'b','LineWidth',1)
% end
% xlabel('Distance upstream from basin outlet, km')
% ylabel('Elevation, m')
% xlim([0 max(Dist./1000)])
% 
% %% Profile of given value
% 
% %SET STARTING LINK
% %loc=1;%001,076
% %loc=42;%074
% %loc=161;%080
% %loc=9;%071
% loc=1;%061
% 
% figure; hold on; box on
% %color=['b','r','k'];
% %yr=[0,1,6];
% for j=1:timesteps-1
% %t=1+365*yr(j);
% t=j;
% 
% index=Connect(loc,~isnan(Connect(loc,:)));
% 
% % STEP THROUGH EACH OF THESE
% y=seddepth(t,index); if j==1;ylabel('Total sediment depth, m^3');end
% %y=seddepthdf(t,index); if j==1;ylabel('Accumulation from debris flow sediment, m^3');end
% %y=Dg(t,index); if j==1;ylabel('Dg (includes sand), m');set(gca,'YScale','log');end
% %y=DG(t,index); if j==1;ylabel('DG (excludes sand), m');set(gca,'YScale','log');end
% %y=Fs(t,index); if j==1;ylabel('Surface sand fraction (Fs), m');end
% 
% %old
% %y=GridID(index,1);
% %[y,yt]=max(seddepthdf(:,index),[],1);
% %y=time(yt);
% %y=Slope(index);
% %y=max(Fs(:,index),[],1);
% %y=min(Dg(:,index),[],1);
% %y(isnan(y))=0;
% %y(seddepth(t,index)<0.1)=NaN;
% 
% x=Dist(index)./1000;%km
% 
% xx=[];
% yy=[];
% % xx=cat(1,xx,0,x(1));
% % yy=cat(1,yy,y(1),y(1));
% % for i=2:length(x)
% %     xx=cat(1,xx,x(i-1:i,1));
% %     yy=cat(1,yy,y(i),y(i));
% % end
% for i=1:length(x)-1
%     xx=cat(1,xx,x(i:i+1,1));
%     yy=cat(1,yy,y(i),y(i));
% end
% xx=cat(1,xx,x(end),0);
% yy=cat(1,yy,y(end),y(end));
% 
% %plot(xx,yy,color(j));
% 
% if j==1
%     plot(xx,yy,'b','LineWidth',2);
% elseif j==1+365
%     plot(xx,yy,'k','LineWidth',2);
% elseif j==timesteps-1
%     plot(xx,yy,'r','LineWidth',2);
% else
%     plot(xx,yy)
% end
% 
% end
% 
% xlabel('Distance upstream from basin outlet, km')
% xlim([0 max(xx)])
% %legend('2 years','10 years','19 years')
% 
% %% Elevation profile -- pathway initial and final
% figure; hold on; box on
% plot([Dist(index)./1000;0],[Elev(1,index),mnelev(index(end),1)],'b','LineWidth',1)
% plot([Dist(index)./1000;0],[Elev(timesteps-1,index),mnelev(index(end),1)],'r','LineWidth',1)
% xlabel('Distance upstream from basin outlet, km')
% ylabel('Elevation, m')
% legend('Initial','Final')
% xlim([0 max(Dist(index)./1000)])
% 
% %% Slope profile -- pathway initial and final
% figure; hold on; box on
% x=Dist(index)./1000;%km
% y1=Slope_raw(index);
% y2=Slope(index);
% xx=[];
% yy1=[];
% yy2=[];
% for i=1:length(x)-1
%     xx=cat(1,xx,x(i:i+1,1));
%     yy1=cat(1,yy1,y1(i),y1(i));
%     yy2=cat(1,yy2,y2(i),y2(i));
% end
% xx=cat(1,xx,x(end),0);
% yy1=cat(1,yy1,y1(end),y1(end));
% yy2=cat(1,yy2,y2(end),y2(end));
% 
% plot(xx,yy1,'b','LineWidth',1);
% plot(xx,yy2,'r','LineWidth',1);
% 
% xlabel('Distance upstream from basin outlet, km')
% ylabel('Bed slope')
% legend('Initial','Final')
% xlim([0 max(Dist(index)./1000)])
% %set(gca,'YScale','log')
% 
% 
% %%
% %% Volume within/exported from network over time
% %%
% %% Volume of debris flow sediment by grainsize in network
% 
% lnkvoldfgs(1:timesteps,1:LinkNum,1:gsclass)=NaN;
% for j=1:gsclass
%     for i=1:LinkNum
%         for tt=1:timesteps
%             %only works when no abrasion! otherwise need to bound P_d
%             lnkvoldfgs(tt,i,j)=sum(P_vol{tt,i}(and(P_idx{tt,i}>=pidxdf,P_d{tt,i}==Dpsd(j))));
%         end
%     end
% end
% 
% %mineral volume, not including porosity
% volinbasin=sum(lnkvoldfgs,2);
% totalvolinbasin=sum(volinbasin,3);
% 
% %% Volume of debris flow sediment by grainsize that has left the network
% %mineral volume, not including porosity
% voloutbasin=repmat(volinbasin(1,1,:),timesteps,1,1)-volinbasin(:,1,:);
% totalvoloutbasin=sum(voloutbasin,3);
% 
% %% Fraction within network
% figure; hold on; box on
% for i=1:gsclass-1
% plot((1:timesteps)./365,volinbasin(:,1,i)./volinbasin(1,1,i))
% end
% ylabel('Fraction of debris flow input within network')
% xlabel('Time since debris flow input, years')
% plot((1:timesteps)./365,totalvolinbasin(:,1,1)./totalvolinbasin(1,1,1),'k','LineWidth',2)
% legend(cat(1,num2str(Dpsd(1:end-1).*1000),' all'))
% ylim([0,1])
% xlim([0 timesteps./365])
% 
% %% Volume within network
% figure; hold on; box on
% for i=1:gsclass-1
% plot((1:timesteps)./365,volinbasin(:,1,i))
% end
% ylabel('Volume of debris flow input within network (no porosity), m^3')
% xlabel('Time since debris flow input, years')
% plot((1:timesteps)./365,totalvolinbasin(:,1,1),'k','LineWidth',2)
% legend(cat(1,num2str(Dpsd(1:end-1).*1000),' all'))
% %ylim([0,1])
% xlim([0 timesteps./365])
% %set(gca,'YScale','log')
% 
% %% Fraction within reservoir
% figure; hold on; box on
% for i=1:gsclass-1
% plot((1:timesteps)./365,voloutbasin(:,1,i)./volinbasin(1,1,i))
% end
% ylabel('Fraction of debris flow input within reservoir')
% xlabel('Time since debris flow input, years')
% plot((1:timesteps)./365,totalvoloutbasin(:,1,1)./totalvolinbasin(1,1,1),'k','LineWidth',2)
% legend(cat(1,num2str(Dpsd(1:end-1).*1000),' all'))
% ylim([0,1])
% xlim([0 timesteps./365])
% 
% %% Volume within reservoir
% figure; hold on; box on
% for i=1:gsclass-1
% plot((1:timesteps)./365,voloutbasin(:,1,i))
% end
% ylabel('Volume of debris flow input within reservoir (no porosity), m^3')
% xlabel('Time since debris flow input, years')
% plot((1:timesteps)./365,totalvoloutbasin(:,1,1),'k','LineWidth',2)
% legend(cat(1,num2str(Dpsd(1:end-1).*1000),' all'))
% %ylim([0,1])
% xlim([0 timesteps./365])
% %set(gca,'YScale','log')
% 
% %% Total volumes
% %total volume DF inputs: network + RES
% totalvolinbasin(1,1,1)%m3
% 
% %total volume DF inputs delivered to RES
% totalvoloutbasin(end,1,1)%m3
% 
% %% Surface area of network
% sum(Length(RES~=1).*B(RES~=1))
% 
% %%
% % figure
% % plot(volinbasin(:,1,7))
% % %%
% % volinbasin(end,1,:)./volinbasin(1,1,:)
% %%
% % lnkvoldfgs_burn=lnkvoldfgs;
% % lnkvoldfgs_burn(:,burn==0,:)=0;
% % volinburn=sum(lnkvoldfgs_burn,2);
% % %%
% % figure; hold on; box on
% % for i=1:gsclass
% % plot((1:length(Qgage))./365,volinburn(:,1,i)./volinburn(1,1,i))
% % end
% % ylabel('Fraction of DF input within burned area')
% % xlabel('Time since debris flow input, years')
% % legend(num2str(Dpsd))
% % %%
% % % volinburn(end,1,:)./volinburn(1,1,:)
% 
% 
% %%
% %% Surface grain size distribution
% %%
% %% Surface grain size distribution
% t=1
% i=1
% %Dft=[2,4,8,32,64,128,256,512]'./1000;%m
% cpsd(1:timesteps,1:gsclass,1)=NaN;
% for t=1:timesteps
% %for i=1:LinkNum
% if isempty(P_loc{t,i})
%     continue
% end
% 
% actidx=P_storage{t,i}==0;
% actvol=sum(P_vol{t,i}(actidx));
% 
% for j=1:gsclass
%     actjidx=and(P_storage{t,i}==0,P_d{t,i}<=Dcsd(j));
%     actftjvol=sum(P_vol{t,i}(actjidx));
%     cpsd(t,j)=actftjvol./actvol;
% end
% end
% %
% figure; hold on; box on
% jj=1:10:(timesteps-2);
% ci=linspace(0.8,0.2,timesteps); 
% for jjj=1:length(jj)
% j=jj(jjj);
% plot(Dcsd,cpsd(j,:),'Color',[ci(j) ci(j) ci(j)])
% end
% plot(Dcsd,cpsd(1+365*0,:),'b','LineWidth',3)
% plot(Dcsd,cpsd(1+365*1,:),'r','LineWidth',3)
% plot(Dcsd,cpsd(1+365*6,:),'k','LineWidth',3)
% set(gca,'XScale','log')
% ylim([0 1])
% 
% % inisspsd=cumsum(Fsspsd);
% % inisfpsd=cumsum(Fsfpsd);
% % inidfpsd=cumsum(Fdfpsd);
% % plot(Dft,inisspsd,'g','LineWidth',2)
% % plot(Dft,inisfpsd,'g','LineWidth',2)
% % plot(Dft,inidfpsd,'g','LineWidth',2)
% title(num2str(i))
% 
% 
% %%
% %% Travel distances and times
% %%
% %% parcel travel distances and travel times
% DF_idx=[];
% DF_d=[];
% DF_loc0=[];
% 
% for i=1:LinkNum
% 
%     dfcellidx=find(P_idx{1,i}>=pidxdf)';
%     if isempty(dfcellidx)
%         continue
%     end
%     DF_idx=cat(1,DF_idx,P_idx{1,i}(dfcellidx)');
%     DF_d=cat(1,DF_d,P_d{1,i}(dfcellidx)');
%     DF_loc0=cat(1,DF_loc0,repmat(i,length(dfcellidx),1));
%     
% end
% %%
% DF_loc7(1:length(DF_idx),1)=NaN;
% DF_tt(1:length(DF_idx),1)=NaN;
% DF_ll(1:length(DF_idx),1)=NaN;
% 
% for dfp=1:length(DF_idx)
% dfp
% index=Connect(DF_loc0(dfp),~isnan(Connect(DF_loc0(dfp),:)));
% 
% t=1;
% i=1;
% while isnan(DF_loc7(dfp))
%     if ~isempty(find(DF_idx(dfp)==P_idx{t,index(i)},1))
%         t=t+1;
%     else
%         i=i+1;
%     end
% 
%     if i>length(index)
%         DF_loc7(dfp)=index(i-1);
%         DF_tt(dfp)=t;      
%     end
%     
%     if t>timesteps
%         DF_loc7(dfp)=index(i);
%         DF_tt(dfp)=t-1;
%     end
%     
%     if i>length(index) && t<timesteps
%         DF_loc7(dfp)=index(i-1);
%         DF_tt(dfp)=t-1;
%     end
% end
% tosum=Connect(DF_loc0(dfp),1:find(Connect(DF_loc0(dfp),:)==DF_loc7(dfp),1));
% if numel(tosum)==1
%     DF_ll(dfp)=Length(tosum(1)).*...
%         P_loc{DF_tt(dfp),DF_loc7(dfp)}(find(P_idx{DF_tt(dfp),DF_loc7(dfp)}==DF_idx(dfp),1));
% elseif isempty(find(P_idx{DF_tt(dfp),DF_loc7(dfp)-1}==DF_idx(dfp),1))
%     DF_ll(dfp)=sum(Length(tosum(1:end)));
% else
%     DF_ll(dfp)=sum(Length(tosum(1:end-1)))+Length(tosum(end)).*...
%         P_loc{DF_tt(dfp),DF_loc7(dfp)}(find(P_idx{DF_tt(dfp),DF_loc7(dfp)-1}==DF_idx(dfp),1));
% end
% end
% %%
% i=5;
% figure
% hist(log10(DF_ll(and(DF_d==Dpsd(i),DF_ll>0))),30)
% xlim([-5 5])
% title(num2str(Dpsd(i)))
% 
% %%
% figure; hold on; box on
% for i=1:gsclass
% x=sort(DF_ll(and(DF_d==Dpsd(i),DF_ll>0)));
% y=(1:length(x))'./length(x);
% plot(x,y)
% end
% set(gca,'XScale','log')
% legend(num2str(Dpsd))
% 
% %%
% DF_tt(DF_tt<timesteps)
% DF_d(DF_tt<timesteps)
% 
% %%
% figure
% plot(DF_ll,DF_tt,'.b')

%%
%% OLD
%%
%%
figure; hold on; box on
plot((1:timesteps)./365,seddepth(:,3),'b')
plot((1:timesteps)./365,seddepth(:,4),'g')
plot((1:timesteps)./365,seddepth(:,10),'r')
%plot((1:length(Q))./365,Q(:,1))
%plot((1:length(Q))./365,seddepth,'b')
%set(gca,'YScale','log')
%ylabel('Total sediment depth within a link, m')
% xlabel('Time since flow input, years')
% xlim([0 70])
% % plot((1:length(Q))./365,seddepth(:,316),'b')
% % plot((1:length(Q))./365,seddepth(:,438),'r')
% % plot((1:length(Q))./365,seddepth(:,416),'k')
% % plot((1:length(Qgage))./365,seddepthdf(:,161),'b')
% % plot((1:length(Qgage))./365,seddepthdf(:,380),'r')
% % plot((1:length(Qgage))./365,seddepthdf(:,424),'k')
% %plot((1:length(Qgage))./365,seddepth(:,267),'k')
% legend('Upstream Link-20','Middle Link-44','Downstream Link-67')
% ylabel('Sediment Depth (m)')
% xlabel('Time since flow input, years')

%%
% figure;
% plot((1:timesteps)./365,seddepth,'b')
% title('most upstream link sed depth (m)')
% %ylim([0 3])
% xlabel('years')
% %%
% figure; hold on; box on
% %plot((1:length(Qgage))./365,Dg)
% %set(gca,'YScale','log')
% ylabel('D50 within a link, m')
% xlabel('Time since flow input, years')
% xlim([0 70])
% plot((1:timesteps)./365,Dg(:,20),'b')
% plot((1:timesteps)./365,Dg(:,38),'r')
% plot((1:timesteps)./365,Dg(:,67),'k')
% legend('Link=20','Link=38','Link=67')
% %%
% figure; hold on; box on
% %plot((1:length(Qgage))./365,Fs)
% %set(gca,'YScale','log')
% ylabel('Sand fraction')
% xlabel('Time since flow input, years')
% xlim([0 70])
% plot((1:timesteps)./365,Fs(:,20),'b')
% plot((1:timesteps)./365,Fs(:,38),'r')
% plot((1:timesteps)./365,Fs(:,67),'k')
% % plot((1:length(Qgage))./365,Fs(:,316),'b')
% % plot((1:length(Qgage))./365,Fs(:,438),'r')
% legend('Link=20','Link=38','Link=67')
% 
% %%
% figure
% plot(P_loc{end,1},'.b')
% figure
% plot(P_loc{1,1},'.r')
% 
% %% contour
% figure; hold on; box on
% loc=76;%middle
% index=Connect(loc,~isnan(Connect(loc,:)));
% y=seddepthdf(:,index);%
% x=cumsum(Length(index))./1000;%km
% contourf(repmat(x',timesteps,1),repmat(time,1,length(x)),(y),'LineStyle','none')
% ylabel('Time since debris flow input, years')
% xlabel('Distance downstream, km')
%%
%%Compute distance from every link to the outlet
Dist(1:LinkNum,1)=NaN;
for i=1:LinkNum
    lastConn=find(Connect(i,:)==OutletLinkID);
    Dist(i,1)=sum(Length(Connect(i,1:lastConn)),1);
end

%% Elevation profile -- all links
figure; hold on; box on
for i=1:LinkNum
    index=Connect(i,~isnan(Connect(i,:)));
    plot([Dist(index,1)./1000; 0],...
        [mxelev(index,1); mnelev(index(end),1)],'b','LineWidth',1)
end
xlabel('Distance upstream from basin outlet, km')
ylabel('Elevation, m')
xlim([0 max(Dist./1000)])

%% Profile of given value

%SET STARTING LINK
%loc=1;%001,076
%loc=42;%074
%loc=161;%080
%loc=9;%071
loc=345;%061

figure; hold on; box on
%color=['b','r','k'];
%yr=[0,1,6];
for j=1:timesteps-1
%t=1+365*yr(j);
t=j;

index=Connect(loc,~isnan(Connect(loc,:)));

% STEP THROUGH EACH OF THESE
y=seddepth(t,index); if j==1;ylabel('Total sediment depth, m');end
%y=seddepthdf(t,index); if j==1;ylabel('Accumulation from debris flow sediment, m');end
%y=Dg(t,index); if j==1;ylabel('Dg (includes sand), m');set(gca,'YScale','log');end
%y=DG(t,index); if j==1;ylabel('DG (excludes sand), m');set(gca,'YScale','log');end
%y=Fs(t,index); if j==1;ylabel('Surface sand fraction (Fs), m');end

%old
%y=GridID(index,1);
%[y,yt]=max(seddepthdf(:,index),[],1);
%y=time(yt);
%y=Slope(index);
%y=max(Fs(:,index),[],1);
%y=min(Dg(:,index),[],1);
%y(isnan(y))=0;
%y(seddepth(t,index)<0.1)=NaN;

x=Dist(index)./1000;%km

xx=[];
yy=[];
% xx=cat(1,xx,0,x(1));
% yy=cat(1,yy,y(1),y(1));
% for i=2:length(x)
%     xx=cat(1,xx,x(i-1:i,1));
%     yy=cat(1,yy,y(i),y(i));
% end
for i=1:length(x)-1
    xx=cat(1,xx,x(i:i+1,1));
    yy=cat(1,yy,y(i),y(i));
end
xx=cat(1,xx,x(end),0);
yy=cat(1,yy,y(end),y(end));

%plot(xx,yy,color(j));

if j==1
    plot(xx,yy,'b','LineWidth',2);
elseif j==1+365/5
    plot(xx,yy,'k','LineWidth',2);
elseif j==timesteps-1
    plot(xx,yy,'r','LineWidth',2);
else
    plot(xx,yy)
end

end

xlabel('Distance upstream from basin outlet, km')
xlim([0 max(xx)])
%legend('2 years','10 years','19 years')