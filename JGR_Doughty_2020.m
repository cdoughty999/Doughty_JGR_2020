%tree girdle experiment
clear
load('JGR_Doughty_2020data.mat');

%% weather Figure 2
datetemp = allweatherdatP1.date(2500:end);
timetemp = allweatherdatP1.time(2500:end);
temp1 = allweatherdatP1.Soiltemperature6cm(2500:end);
temp2 = allweatherdatP2.Soiltemperature6cm(2500:end);
temp3 = allweatherdatP3.Soiltemperature6cm(2500:end);


temp1air = allweatherdatP1.Airtemperature12cm(2500:end);
temp2air = allweatherdatP2.Airtemperature12cm(2500:end);
temp3air = allweatherdatP3.Airtemperature12cm(2500:end);

datesm = allweatherdatP1.date(2500:end);
timesm = allweatherdatP1.time(2500:end);
sm1 = allweatherdatP1.Volmoisture(2500:end);
sm2 = allweatherdatP2.Volmoisture(2500:end);
sm3 = allweatherdatP3.Volmoisture(2500:end);
% subplot(211)
% plot(sm1-sm3)
% subplot(212)
% plot(temp1-temp3)
t=1;
for i= 0:.04:.96
    ix = find(timetemp>i & timetemp<(i+.05));
    temp1a(t)= nanmean(temp1(ix));
    temp2a(t)= nanmean(temp2(ix));
    temp3a(t)= nanmean(temp3(ix));
    temp1as(t)= nanstd(temp1(ix));
    temp2as(t)= nanstd(temp2(ix));
    temp3as(t)= nanstd(temp3(ix));
    
    temp1aira(t)= nanmean(temp1air(ix));
    temp2aira(t)= nanmean(temp2air(ix));
    temp3aira(t)= nanmean(temp3air(ix));
    temp1airas(t)= nanstd(temp1air(ix));
    temp2airas(t)= nanstd(temp2air(ix));
    temp3airas(t)= nanstd(temp3air(ix));
    t=t+1;
end
v=6;
temp1az = [temp1a(v:end) temp1a(1:v-1)];
temp2az = [temp2a(v:end) temp2a(1:v-1)];
temp3az = [temp3a(v:end) temp3a(1:v-1)];
temp1asz = [temp1as(v:end) temp1as(1:v-1)];
temp2asz = [temp2as(v:end) temp2as(1:v-1)];
temp3asz = [temp3as(v:end) temp3as(1:v-1)];

temp1airaz = [temp1aira(v:end) temp1aira(1:v-1)];
temp2airaz = [temp2aira(v:end) temp2aira(1:v-1)];
temp3airaz = [temp3aira(v:end) temp3aira(1:v-1)];
temp1airasz = [temp1airas(v:end) temp1airas(1:v-1)];
temp2airasz = [temp2airas(v:end) temp2airas(1:v-1)];
temp3airasz = [temp3airas(v:end) temp3airas(1:v-1)];

tt=1;
for i= 0:.04:.96
    ix = find(timesm>i & timesm<(i+.05));
    sm1a(tt)= nanmean(sm1(ix));
    sm2a(tt)= nanmean(sm2(ix));
    sm3a(tt)= nanmean(sm3(ix));
    sm1as(tt)= nanstd(sm1(ix));
    sm2as(tt)= nanstd(sm2(ix));
    sm3as(tt)= nanstd(sm3(ix));
    tt=tt+1;
end
sm1az = [sm1a(v:end) sm1a(1:v-1)];
sm2az = [sm2a(v:end) sm2a(1:v-1)];
sm3az = [sm3a(v:end) sm3a(1:v-1)];
sm1asz = [sm1as(v:end) sm1as(1:v-1)];
sm2asz = [sm2as(v:end) sm2as(1:v-1)];
sm3asz = [sm3as(v:end) sm3as(1:v-1)];


figure(1)
subplot(511)
errorbar(temp1az, temp1asz/sqrt(24), 'Color', [17 17 17],'LineWidth',2)
hold on
%errorbar(temp2az, temp2asz/sqrt(24),'r')
errorbar(temp3az, temp3asz/sqrt(24),'k','LineWidth',2)
xlabel('time')
ylabel('Soil Temp (C)')
legend('Plot T-17', 'Plot NT')
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(1),ylim(2),'A')
subplot(512)
errorbar(temp1airaz, temp1airasz/sqrt(24),'Color', [17 17 17],'LineWidth',2)
hold on
errorbar(temp3airaz, temp3airasz/sqrt(24),'k','LineWidth',2)
xlabel('time')
ylabel('Air Temp (C)')
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(1),ylim(2),'B')
subplot(513)
errorbar(sm1az, sm1asz/sqrt(24),'LineWidth',3)
hold on
% errorbar(sm2az, sm2asz/sqrt(24),'r')
errorbar(sm3az, sm3asz/sqrt(24),'g','LineWidth',3)
xlabel('time')
ylabel('Volumetric water %')
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(1),ylim(2),'C')
subplot(514)
plot(temp1air)
hold on
% errorbar(sm2az, sm2asz/sqrt(24),'r')
plot(temp3air,'g')
set(gca,'xtick',[])
%gtext(['Oct 5, 2019']);
%gtext(['Aug 26, 2020']);

ylabel('Air Temp at 12cm')
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(1),ylim(2),'D')
subplot(515)
plot(sm1)
hold on
% errorbar(sm2az, sm2asz/sqrt(24),'r')
plot(sm3,'g')
set(gca,'xtick',[])
%gtext(['Oct 5, 2019']);
%gtext(['Aug 26, 2020']);

ylabel('Volumetric water %')
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(1),ylim(2),'E')
%% figure 3 NPP_all
%NPP_AWC, followed by NPP_litter and NPP_fine roots. 
figure(2)
subplot(411)
errorbar(NPP_ACW(:,1), NPP_ACW(:,2),'r','LineWidth',3);
hold on
errorbar(NPP_ACW(:,3), NPP_ACW(:,4),'g','LineWidth',3);
errorbar(NPP_ACW(:,5), NPP_ACW(:,6),'LineWidth',3);
xlabel('Month')
ylabel('Mg C Ha^-^1 mo^-^1')
legend('Plot T-17', 'Plot T-15', 'Plot NT')
title('NPP ACW')
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(1),ylim(2),'A')
subplot(412)
errorbar(NPP_litterfall(:,1), NPP_litterfall(:,2)./sqrt(6),'r','LineWidth',3);
hold on
errorbar(NPP_litterfall(:,3), NPP_litterfall(:,4)./sqrt(6),'g','LineWidth',3);
errorbar(NPP_litterfall(:,5), NPP_litterfall(:,6)./sqrt(6),'LineWidth',3);
xlabel('Month')
ylabel('Mg C Ha^-^1 mo^-^1')
title('NPP litterfall')
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(1),ylim(2),'B')
subplot(413)
errorbar(NPP_fineroot(:,1), NPP_fineroot(:,2)./sqrt(3),'r','LineWidth',3);
hold on
errorbar(NPP_fineroot(:,3), NPP_fineroot(:,4)./sqrt(3),'g','LineWidth',3);
errorbar(NPP_fineroot(:,5), NPP_fineroot(:,6)./sqrt(3),'LineWidth',3);
xlabel('Month')
ylabel('Mg C Ha^-^1 mo^-^1')
title('NPP fine root')
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(1),ylim(2),'C')
subplot(414)
errorbar(NPP_all(:,1), NPP_all(:,2)./sqrt(6),'r','LineWidth',3);
hold on
errorbar(NPP_all(:,3), NPP_all(:,4)./sqrt(6),'g','LineWidth',3);
errorbar(NPP_all(:,5), NPP_all(:,6)./sqrt(6),'LineWidth',3);
xlabel('Month')
ylabel('Mg C Ha^-^1 mo^-^1')
title('NPPtotal')
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(1),ylim(2),'D')
%% total soil respiration
figure(3)
subplot(211)
errorbar(Rsoiltot(:,1), Rsoiltot(:,2)./sqrt(6),'r','LineWidth',3);
hold on
errorbar(Rsoiltot(:,3), Rsoiltot(:,4)./sqrt(6),'g','LineWidth',3);
errorbar(Rsoiltot(:,5), Rsoiltot(:,6)./sqrt(6),'LineWidth',3);
xlabel('Month')
ylabel('Mg C Ha^-^1 mo^-^1')
title('Total Soil Respiration')
legend('Plot T-17', 'Plot T-15', 'Plot NT')
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(1),ylim(2),'A')
subplot(212)
errorbar(RsoilA(:,1), RsoilA(:,2)./sqrt(6),'r','LineWidth',3);
hold on
errorbar(RsoilA(:,3), RsoilA(:,4)./sqrt(6),'g','LineWidth',3);
errorbar(RsoilA(:,5), RsoilA(:,6)./sqrt(6),'LineWidth',3);
xlabel('Month')
ylabel('Mg C Ha^-^1 mo^-^1')
title('R Rhizosphere')
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(1),ylim(2),'B')
%% tree respiration

figure(4)
subplot(211)
errorbar(treeres(:,1), treeres(:,2)./sqrt(6),'r','LineWidth',3);
hold on
errorbar(treeres(:,3), treeres(:,4)./sqrt(6),'g','LineWidth',3);
errorbar(treeres(:,5), treeres(:,6)./sqrt(6),'LineWidth',3);
xlabel('Month')
ylabel('g CO^2 m2 hr^-^1')
title('Mean Tree Wood Respiration')
legend('Plot T-17', 'Plot T-15', 'Plot NT')
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(1),ylim(2),'A')
%convert tree res to MgC ha mo
P1SA = 1847;% m2 wood area from WoodSA Plot 1 G34
P2SA = 1426;% m2 wood area from WoodSA Plot 2 G34
P3SA = 5426;% m2 wood area from WoodSA Plot 3 G34
tempcorrwr = 0.73;% temp correction from Fig5woodres temp
unitcon = (((365*24)/1000000)*(12/28))/12; % convert g CO2 m2 hr-1 to Mg Ha mo 

Rwoodplotp1 = treeres(:,1:2).*(P1SA*tempcorrwr*unitcon);
Rwoodplotp2 = treeres(:,3:4).*(P2SA*tempcorrwr*unitcon);
Rwoodplotp3 = treeres(:,5:6).*(P3SA*tempcorrwr*unitcon);
Rwoodall = [Rwoodplotp1 Rwoodplotp2 Rwoodplotp3];
subplot(212)
errorbar(Rwoodplotp1(:,1), Rwoodplotp1(:,2)./sqrt(6),'r','LineWidth',3);
hold on
errorbar(Rwoodplotp2(:,1), Rwoodplotp2(:,2)./sqrt(6),'g','LineWidth',3);
errorbar(Rwoodplotp3(:,1), Rwoodplotp3(:,2)./sqrt(6),'LineWidth',3);
xlabel('Month')
ylabel('Mg C Ha^-^1 mo^-^1')
title('R Stems')
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(1),ylim(2),'B')
%% Ra, PCE, CUE
%Mean foliage respiration normalized to 10 °C was 0.20 μmol m–2 (hemi‐leaf surface area) s–1
% conversion 0.757 = LAI*((10000*60*60*24*365*12)/1000000000000)
Rcanopyf = Rcanopy.*0.67;% To account for light inhibition of dark respiration we multiplied our result by 0.67(Malhi et al., 2009)
Ra = Rwoodall + RsoilA + Rcanopyf;
PCE = Ra + NPP_all;
CUE = NPP_all./PCE;

figure(5)
subplot(311)
errorbar(Ra(:,1), Ra(:,2)./sqrt(6),'r','LineWidth',3);
hold on
errorbar(Ra(:,3), Ra(:,4)./sqrt(6),'g','LineWidth',3);
errorbar(Ra(:,5), Ra(:,6)./sqrt(6),'LineWidth',3);
xlabel('Month')
ylabel('Mg C Ha^-^1 mo^-^1')
title('Ra')
legend('Plot T-17', 'Plot T-15', 'Plot NT')
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(1),ylim(2),'A')
subplot(312)
errorbar(PCE(:,1), PCE(:,2)./sqrt(6),'r','LineWidth',3);
hold on
errorbar(PCE(:,3), PCE(:,4)./sqrt(6),'g','LineWidth',3);
errorbar(PCE(:,5), PCE(:,6)./sqrt(6),'LineWidth',3);
xlabel('Month')
ylabel('Mg C Ha^-^1 mo^-^1')
title('Plant Carbon Expenditure')
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(1),ylim(2),'B')
subplot(313)
errorbar(CUE(:,1), CUE(:,2)./sqrt(6),'r','LineWidth',3);
hold on
errorbar(CUE(:,3), CUE(:,4)./sqrt(6),'g','LineWidth',3);
errorbar(CUE(:,5), CUE(:,6)./sqrt(6),'LineWidth',3);
xlabel('Month')
ylabel('Mg C Ha^-^1 mo^-^1')
title('Carbon Use Efficiency')
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(1),ylim(2),'C')
%% herbivory

figure(6)
subplot(211)
errorbar(totherb(:,1), totherb(:,2),'r','LineWidth',3);
hold on
errorbar(totherb(:,3), totherb(:,4),'g','LineWidth',3);
errorbar(totherb(:,5), totherb(:,6),'LineWidth',3);
xlabel('Month')
ylabel('Mg C Ha^-^1 mo^-^1')
title('Total Herbivory')
legend('Plot T-17', 'Plot T-15', 'Plot NT')
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(1),ylim(2),'A')
subplot(212)
errorbar(indanimalherb(:,1), indanimalherb(:,2),'r','LineWidth',3);
hold on
errorbar(indanimalherb(:,3), indanimalherb(:,4),'g','LineWidth',3);
errorbar(indanimalherb(:,5), indanimalherb(:,6),'k','LineWidth',3);
errorbar(indanimalherb(:,7), indanimalherb(:,8),'LineWidth',3);
xlabel('Month')
ylabel('Mg C Ha^-^1 mo^-^1')
title('Herbivory by species')
legend('Proghorn', 'Deer', 'Elk', 'Total')
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(1),ylim(2),'B')

%% Table 4
% order NPP_AWC, followed by NPP_litter and NPP_fine rootsNPP_litterfall
%Herbdung is from Herb_est cell U9
%NPP understory from Biomass understory
TotalNPPan = sum(NPP_ACW)+ sum(NPP_litterfall)+ sum(NPP_fineroot)+ sum(NPP_ACW)*0.25+ NPPunder;

%sum errors with quadrature
ACWstd1 = sum(NPP_ACW(:,2));
ACWstd2 = sum(NPP_ACW(:,4));
ACWstd3 = sum(NPP_ACW(:,6));
FRstd1 = sum(NPP_fineroot(:,2));
FRstd2 = sum(NPP_fineroot(:,4));
FRstd3 = sum(NPP_fineroot(:,6));
LFstd1 = sum(NPP_litterfall(:,2));
LFstd2 = sum(NPP_litterfall(:,4));
LFstd3 = sum(NPP_litterfall(:,6));
totNPPstd1 = sqrt(ACWstd1^2 + FRstd1^2 + (ACWstd1*0.25)^2 + LFstd1^2 );
totNPPstd2 = sqrt(ACWstd2^2 + FRstd2^2 + (ACWstd2*0.25)^2 + LFstd2^2 );
totNPPstd3 = sqrt(ACWstd3^2 + FRstd3^2 + (ACWstd3*0.25)^2 + LFstd3^2 );

rsstd1 = sum(RsoilA(:,2));
rsstd2 = sum(RsoilA(:,4));
rsstd3 = sum(RsoilA(:,6));
rwstd1 = sum(Rwoodall(:,2));
rwstd2 = sum(Rwoodall(:,4));
rwstd3 = sum(Rwoodall(:,6));
rcstd1 = sum(Rcanopy(:,2));
rcstd2 = sum(Rcanopy(:,4));
rcstd3 = sum(Rcanopy(:,6));
totrastd1 = sqrt(rsstd1^2 + rwstd1^2 +  rcstd1^2 );
totrastd2 = sqrt(rsstd2^2 + rwstd2^2 +  rcstd2^2 );
totrastd3 = sqrt(rsstd3^2 + rwstd3^2 +  rcstd3^2 );

GPPFstdte1 = sqrt(totNPPstd1^2 + totrastd1^2);
GPPFstdte2 = sqrt(totNPPstd2^2 + totrastd2^2);
GPPFstdte3 = sqrt(totNPPstd3^2 + totrastd3^2);
CUEAstd1 = abs(mean(CUE(:,1))) * sqrt( ((totNPPstd1/sum(NPP_all(:,1)))^2) + ((GPPFstdte1/sum(PCE(:,1)))^2) ); 
CUEAstd2 = abs(mean(CUE(:,3))) * sqrt( ((totNPPstd2/sum(NPP_all(:,3)))^2) + ((GPPFstdte2/sum(PCE(:,3)))^2) ); 
CUEAstd3 = abs(mean(CUE(:,5))) * sqrt( ((totNPPstd3/sum(NPP_all(:,5)))^2) + ((GPPFstdte3/sum(PCE(:,5)))^2) ); 

TotalNPPan2 = [TotalNPPan(1)  totNPPstd1 TotalNPPan(3)  totNPPstd2 TotalNPPan(5) totNPPstd3];
Ra2 = [sum(Ra(:,1)) totrastd1 sum(Ra(:,3)) totrastd1 sum(Ra(:,5)) totrastd1];
PCE2 = [sum(PCE(:,1)) GPPFstdte1 sum(PCE(:,3)) GPPFstdte2 sum(PCE(:,5)) GPPFstdte3];
CUE2 = [mean(CUE(:,1)) CUEAstd1 mean(CUE(:,3)) CUEAstd2 mean(CUE(:,5)) CUEAstd3];

Table4 = [sum(NPP_ACW); sum(NPP_litterfall); sum(NPP_fineroot); sum(NPP_ACW)*0.25; NPPunder; TotalNPPan2; sum(totherb); Herbdung; sum(Rwoodall); sum(Rcanopyf); sum(RsoilA); Ra2; sum(Rsoilhet); PCE2; CUE2];


%%Table 5

allwood = [sum(NPP_ACW(:,1))/TotalNPPan(1) sum(NPP_ACW(:,3))/TotalNPPan(3) sum(NPP_ACW(:,5))/TotalNPPan(5)];
alllitter = [sum(NPP_litterfall(:,1))/TotalNPPan(1) sum(NPP_litterfall(:,3))/TotalNPPan(3) sum(NPP_litterfall(:,5))/TotalNPPan(5)];
allroot = [sum(NPP_fineroot(:,1))/TotalNPPan(1) sum(NPP_fineroot(:,3))/TotalNPPan(3) sum(NPP_fineroot(:,5))/TotalNPPan(5)];

Rallwood = [sum(Rwoodall(:,1))/Ra2(1) sum(Rwoodall(:,3))/Ra2(3) sum(Rwoodall(:,5))/Ra2(5)];
Ralllitter = [sum(Rcanopyf(:,1))/Ra2(1) sum(Rcanopyf(:,3))/Ra2(3) sum(Rcanopyf(:,5))/Ra2(5)];
Rallroot = [sum(RsoilA(:,1))/Ra2(1) sum(RsoilA(:,3))/Ra2(3) sum(RsoilA(:,5))/Ra2(5)];

herballdung = [Herbdung(1)./NPPunder(1) Herbdung(3)./NPPunder(3) Herbdung(5)./NPPunder(5)];
herballcam = [sum(totherb(:,1))./NPPunder(1) sum(totherb(:,3))./NPPunder(3)  sum(totherb(:,5))./NPPunder(5)];
GPPherbdung = [Herbdung(1)./PCE2(1) Herbdung(3)./PCE2(3) Herbdung(5)./PCE2(5)];
%from Fig7_cameratraps cell H36,G96, G54
carnherb = [.0004 0.0003 0.0006];

Table5 = [allwood; alllitter; allroot; Rallwood; Ralllitter; Rallroot; herballdung; herballcam; carnherb; GPPherbdung].*100;

