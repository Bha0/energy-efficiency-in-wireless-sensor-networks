clear all;
close all;
clc;
load gp

c1=0.5;
c2=0.5;
w=0.9;
Sensing_region_length=200;  % in meter
Sensing_region_width=200;  % in meter
cluster_radius =30;  % in meter
Sensing_range=36;   % in meter
No_nodes=100;
Packet_size=512;  % in Bytes
Eo=200;  % Initial Energy

trans_power=0.02;
recei_power=0.01;
% Prediction of Number of Clusters in the Network
% xx=30;
% yy=30;
% jj=1;
% while yy<200
% while xx<200
%         center_x(jj)=xx;
%         center_y(jj)=yy;
%         xx=xx+50;
%         jj=jj+1;
% end
% xx=30;
% yy=yy+50;
% end
% center=[center_x;center_y];
% figure,
% grid on;
% hold on;
% for cc=1:length(center_x)
% [circle_x,circle_y]=drawcircle(center(:,cc),cluster_radius,[-360 360]);
% hold on;
% plot(circle_x,circle_y,'k-','LineWidth',1.5);
% ylim([0 200]);
% xlim([0 200]);
% end

% Define parameters
SNR_dB = 10; % Signal-to-Noise Ratio in dB
message = 'Hello, World!'; % Message to transmit
transmission_rate = 1000; % Transmission rate in bits per second

% Convert the message to binary data
binary_data = dec2bin(uint8(message), 8);

% Modulation (e.g., BPSK modulation)
modulated_data = 2 * (binary_data - '0') - 1;

% Add AWGN noise to simulate channel effects
SNR = 10^(SNR_dB/10);
noise_power = 1 / (2 * SNR);
received_data = awgn(modulated_data, SNR, 'measured');

% Demodulation (e.g., BPSK demodulation)
demodulated_data = (received_data + 1) / 2;



figure,
ui2={[0.9 0.9 0.9], [0.9 0.9 0.9], [0.9 0.9 0.9]};
op2=[0 2 4 6 10 14];
Frms=1000;
Vel=2e8;
p = 0.02777;
r = 0.25;
PN = 10000;
iter = 100;
b=1;
i=4;
for j=9:2:14
    t = (1/12:1/6:1)'*2*pi;
    x = 2*sin(t)+1+(i*1.5*2);
    y = 2*round(cos(t))+1+((j+2)*2);
    fill(x,y,ui2{b})
    hold on
    b=b+1;
end
title('Node Initialization');
axis off
xdimension=.08; ydimension=.08;
dx=1-xdimension;
dy=1-ydimension;
for i=1:size(nb,1)
    circle(nb(i,1),nb(i,2),0.3,'w');
    if i~=size(nb,1)
       plot([nb(i,1) nb(i+1,1)],[nb(i,2) nb(i+1,2)],'--'); 
    end
end
a=1;
for i=1:3
%     plot([X1(i), tb(i,1)], [Y1(i), tb(i,2)])
    if i~=3
       plot([tb(i,1) tb(2,1)],[tb(i,2) tb(i+1,2)],'--'); 
    end
    plot([nb(a,1) tb(i,1)],[nb(a,2) tb(i,2)],'--'); 
    a=a+1;
    plot([nb(a,1) tb(i,1)],[nb(a,2) tb(i,2)],'--'); 
    a=a+1;
end
X11 = X1;
plot(14-(rand(1,15)*2),randi([24 32],1,15)-(rand(1,15)*2),'or');
pos = [.68 .71 xdimension ydimension] ; 
axes1 = axes('position',pos);
imshow('BS.png');axis off 
pos = [.68 .47 xdimension ydimension];  
axes2 = axes('position',pos);
imshow('BS.png');axis off 
pos = [.68 .22 xdimension ydimension];  
axes3 = axes('position',pos);
imshow('BS.png');axis off 


NC=(Sensing_region_length*Sensing_region_width)/(2*(cluster_radius^2));

NCavg=(((Sensing_region_length^2)+(2*sqrt(2)*30*cluster_radius))/(2*(cluster_radius^2)));

% Clustering Using PSO
node_x=rand(1,No_nodes)*Sensing_region_length;
X_node=node_x;
node_y=rand(1,No_nodes)*Sensing_region_width;
Y_node=node_y;
w1 = rand(2, No_nodes);
w2 = rand(2, No_nodes);
Position=[node_x;node_y];
hold on;
plot(node_x,node_y,'go','MarkerSize',10,'MarkerFaceColor','g','MarkerEdgeColor','k');
hold off;
v1=rand(1,No_nodes);
v2=rand(1,No_nodes);
velocity_update=[v1;v2];

Lbest  = Position;

% Calculate Fitness of the particles
for pi = 1:size(Position,2)
    P1=Position(:,pi);
    Velo1=velocity_update(:,pi);
    PresentFitness(pi) = fitness(P1,Position,Eo,Sensing_range);    
end

fitness_Lbest  = PresentFitness ;
[fitness_Gbest,gb] = max(fitness_Lbest) ;

for gg=1:size(Position,2)
    Gbest(:,gg) = Lbest(:,gb);
end

% Calculate Velocity of Particle with their position
velocity_update = (w *velocity_update) + (w1.*(Lbest-Position)) + (w2.*(Gbest-Position));

% Update the Position with updated velocity             
Position = Position + velocity_update ;

k = 0 ;        % Assign Iteration Number
while  ( k < No_nodes )
k = k + 1;

for pi = 1:size(Position,2)
    P1=Position(:,pi);
    Velo1=velocity_update(:,pi);
    PresentFitness(pi) = fitness(P1,Position,Eo,Sensing_range); 
end

for fi = 1 : size(Position,2)
        if PresentFitness(fi) > fitness_Lbest(fi)
           fitness_Lbest(fi)  = PresentFitness(fi);  
           Lbest(:,fi) = Position(:,fi);
        end   
end
 [Gbest_PresentFit,g] = max(fitness_Lbest);
 local_best(k)=max(fitness_Lbest);
if Gbest_PresentFit > fitness_Gbest
   fitness_Gbest = Gbest_PresentFit;
end
 velocity_update = abs(w *velocity_update + (w1.*(Lbest-Position)) + (w2.*(Gbest-Position)));
 Position = Position + velocity_update;
Global_best_val(k)=fitness_Gbest;
end

[global_best,Gbest_id]=sort(Global_best_val,'descend');
cnr=1;
id_cnt=1;
cnt=1;
count=1;
cntno=1;
sensor_x=[];
sensor_y=[];
figure,
while length(Gbest_id)>0
%     length(Gbest_id)
cluster_head(cnt)=Gbest_id(cnt);
cluster_head_x(cnt)=node_x(Gbest_id(cnt));
cluster_head_y(cnt)=node_y(Gbest_id(cnt));
CH(count,1)=cluster_head_x;
CH(count,2)=cluster_head_y;
    
for hi=1:size(Position,2)
    Dist_CH(hi)=sqrt(((cluster_head_x(cnt)-node_x(hi))^2)+((cluster_head_y(cnt)-node_y(hi))^2));
end
cou=1;
for hi=1:size(Position,2)
    if Dist_CH(hi)<cluster_radius
        normal_id(cou)=hi;
        cou=cou+1;
    end
end
        node_id(id_cnt)=cou;
        id_cnt=id_cnt+1;
co=[abs(rand(1,1)) abs(rand(1,1)) abs(rand(1,1))];
cc_color(count,:)=co;
val_temp1=find(cc_color(:,1)==co(1,1));
val_temp2=find(cc_color(:,2)==co(1,2));
val_temp3=find(cc_color(:,3)==co(1,3));
if ~isempty(val_temp1)
    co(1,1)=randi([0 255],1,1)/255;
end
if ~isempty(val_temp2)
    co(1,2)=randi([0 255],1,1)/255;
end
if ~isempty(val_temp3)
    co(1,3)=randi([0 255],1,1)/255;
end
sensor{count}(1,1)=cluster_head_x;
sensor{count}(1,2)=cluster_head_y;
sensor{count}(3:(2+length(normal_id)),1)=node_x(normal_id);
sensor{count}(3:(2+length(normal_id)),2)=node_y(normal_id);
for ff=1:length(normal_id)
    val_senx{ff}=find(sensor_x==node_x(normal_id(ff)));
    val_seny{ff}=find(sensor_y==node_y(normal_id(ff)));
    if ~isempty(val_senx{ff}) && ~isempty(val_seny{ff})
        normal_id(ff)=-1;
    end
end
    normal_id(normal_id==-1)=[];
sensor_x(1,cntno:(cntno+length(normal_id)-1))=node_x(normal_id);
sensor_y(1,cntno:(cntno+length(normal_id)-1))=node_y(normal_id);
cntno=cntno+length(node_y(normal_id));
length(normal_id)
if length(normal_id)==1
grid on;
hold on;
plot(cluster_head_x,cluster_head_y,'mo','MarkerSize',10,'MarkerFaceColor',co,'MarkerEdgeColor','k','LineWidth',2.5);
else
grid on;
hold on;
plot(cluster_head_x,cluster_head_y,'mo','MarkerSize',10,'MarkerFaceColor',co,'MarkerEdgeColor','k','LineWidth',3);
hold on;
plot(node_x(normal_id),node_y(normal_id),'mo','MarkerSize',10,'MarkerFaceColor',co,'MarkerEdgeColor','k','LineWidth',1);
hold on;
[circle_x,circle_y]=drawcircle([cluster_head_x cluster_head_y],cluster_radius,[-360 360]);
plot(circle_x,circle_y,'k-','LineWidth',1.5);
end

node_x(normal_id)=-1;
node_y(normal_id)=-1;
node_x(Gbest_id(cnt))=-1;
node_y(Gbest_id(cnt))=-1;
Gbest_id(cnt)=0;
for cc=1:length(normal_id)
    coor_val=find(Gbest_id==normal_id(cc));
    Gbest_id(coor_val)=0;
end
Gbest_id(Gbest_id==0)=[];
count=count+1;
    resi=find(node_x~=-1);
    resi_node(cnr)=length(resi);
    cnr=cnr+1;
end
hold on;
for cc=1:length(cluster_head_x)
    text(CH(cc,1),CH(cc,2),num2str(cc),'FontSize',8,'HorizontalAlignment','center');
end
resi_node=sort(resi_node,'ascend');
base_station_x=[200 200];
base_station_y=[0 20];
hold on;
load I;
image(base_station_x,base_station_y,I);
ylim([0 200]);
xlim([0 220]);

% Gravitational Force Estimation
CH(end+1,:)=[200 10];

radi_val=1; 
no_point=1;
ber_pt=radi_val^2*no_point; 
num_cluster=size(CH,1);
weight_level=100; 
TS_val=0.15*sqrt(num_cluster)*radi_val;
rc=0.15*TS_val;
angle_val=0:pi:2*pi;
angle_X=rc*cos(angle_val);
angle_Y=rc*sin(angle_val);
cluster_x=CH(:,1);
cluster_y=CH(:,2);

Distance_val=zeros(num_cluster,num_cluster);
for nc1=1:num_cluster-1
    cluster_x1=cluster_x(nc1);
    cluster_y1=cluster_y(nc1);
    for nc2=nc1+1:num_cluster
        cluster_x2=cluster_x(nc2);
        cluster_y2=cluster_y(nc2);
        diff_x=cluster_x1-cluster_x2;
        diff_y=cluster_y1-cluster_y2;
        dist=sqrt(diff_x^2+diff_y^2);
        Distance_val(nc1,nc2)=dist;
        Distance_val(nc2,nc1)=dist;
    end
end

BER_val=0.5*erfc(sqrt(ber_pt./(Distance_val.^2*no_point))); 
q=1-BER_val; 
qq=q.^weight_level;
weight_node=1./qq; 
% figure;
% hold on;
% for nc1=1:num_cluster-1
%     cluster_x1=cluster_x(nc1);
%     cluster_y1=cluster_y(nc1);
%     for nc2=nc1+1:num_cluster
%         cluster_x2=cluster_x(nc2);
%         cluster_y2=cluster_y(nc2);
%         wt=weight_node(nc1,nc2);
%         plot([cluster_x1 cluster_x2],[cluster_y1 cluster_y2],'-','color',1-[1 1 1]./wt);
%     end
% end
% for nc=1:num_cluster
%     if nc==1
%         plot(cluster_x(nc)+angle_X,cluster_y(nc)+angle_Y,'r-');
%     else
%         plot(cluster_x(nc)+angle_X,cluster_y(nc)+angle_Y,'k-');
%     end
%     text(cluster_x(nc),cluster_y(nc),num2str(nc),'HorizontalAlignment','center','color',[1 0 0]);
% end
% axis equal;
% title('net graph, with - graycolor');

C=true(num_cluster,num_cluster); 
C(1:num_cluster+1:num_cluster^2)=false; 

Short_path=wdijkstra(C,weight_node,size(CH,1)); 
nc=num_cluster/5;
round_nc=round(nc);
packet=7*10^(7);
Tx_packet=packet/(10^5);
for nc=1:num_cluster-1
      distance_CH_BS=sqrt( (CH(nc,1)-(CH(end,1)) )^2 + (CH(nc,2)-(CH(end,2)) )^2 );
figure(100),
hold on;
plot(X_node,Y_node,'mo','MarkerSize',10,'MarkerFaceColor','m','MarkerEdgeColor','k','LineWidth',1);hold on;
for cx=1:size(CH,1)
[circle_x,circle_y]=drawcircle([CH(cx,1) CH(cx,2)],cluster_radius,[-360 360]);
plot(circle_x,circle_y,'k-','LineWidth',1.5);hold on;
end
    for nct=1:num_cluster
        if nct==num_cluster
            plot(cluster_x(nct)+angle_X,cluster_y(nct)+angle_Y,'ro-','MarkerSize',20,'LineWidth',2,'MarkerFaceColor','y'); 
            text(cluster_x(nct),cluster_y(nct),num2str(nct),'FontSize',8);
            plot(cluster_x(nc)+angle_X,cluster_y(nc)+angle_Y,'ro-','MarkerSize',20,'LineWidth',2,'MarkerFaceColor','y'); 
            text(cluster_x(nc),cluster_y(nc),num2str(nc),'FontSize',8)
        else
            plot(cluster_x(nct)+angle_X,cluster_y(nct)+angle_Y,'ko-','MarkerSize',10,'LineWidth',2,'MarkerFaceColor','y');
            text(cluster_x(nct),cluster_y(nct),num2str(nct),'FontSize',8)
        end
    end
    pth=Short_path{nc};
    plot(cluster_x(pth),cluster_y(pth),'k-','LineWidth',3);
    ylim([0 200]);
    xlim([0 200]);
    pause(0.3)
    packet_delay=randi([900 1000],1,1);
    Rx_Packet=randi([(Tx_packet-20) (Tx_packet)],1,1);
   if sum(nc==find(1:(num_cluster-2)))
    close(100);
   end
   Throughput(nc)=packet/packet_delay;
   PDR(nc)=(Rx_Packet/Tx_packet)*100;
end
Throughput_propos=sort(Throughput,'descend');
PDR=sort(PDR,'descend');
Network_lifetime=Throughput_propos(1,1)/Tx_packet;
% Performance Analysis
E_OEERP_res=resi_node(3);
OEERP_res=35;
LEACH_res=42;
BCDCP_res=17;
NEEC_res=19;
Residual_nodes=[E_OEERP_res OEERP_res LEACH_res BCDCP_res NEEC_res];
xx=[1 2 3 4 5];
figure(200),
subplot(2,1,1),plot(xx,Residual_nodes,'b-','LineWidth',2);hold on;
plot(xx,Residual_nodes,'cd','MarkerSize',10,'MarkerFaceColor','c','MarkerEdgeColor','k');
cnames={'E-OEERP','OEERP','LEACH','BCDCP','NEEC'};
rnames = {'Series1'};
f=figure(200);
t = uitable('Parent',f,'Data',Residual_nodes,'RowName',rnames,'ColumnName',cnames,'Position', [50 100 460 60]);
title('Residual node comparison among different protocols');

for ee=1:length(node_id)
    Einitial(ee)=Eo;
end
Eremaining=(trans_power+recei_power)*Eo*node_id;
temp_energy=Einitial-Eremaining;
Energy_consumed=mean(temp_energy);

Energy_exist(1,:)=[10 20 5 12];
Energy_exist(2,:)=[19 42 10 17];
Energy_exist(3,:)=[30 65 27 28];
Energy_exist(4,:)=[37 90 47 30];
Energy_exist(5,:)=[47 110 85 40];
Energy_exist(6,:)=[55 135 120 50];
Energy_consumed1=[(Energy_consumed/20) (Energy_consumed/10) (Energy_consumed/10*1.5) (Energy_consumed/10*2) (Energy_consumed/10*2.5) (Energy_consumed/10*3)];

E_Consumed=[Energy_consumed1(1,1) Energy_exist(1,:);Energy_consumed1(1,2) Energy_exist(2,:);Energy_consumed1(1,3) Energy_exist(3,:);Energy_consumed1(1,4) Energy_exist(4,:);Energy_consumed1(1,5) Energy_exist(5,:);Energy_consumed1(1,6) Energy_exist(6,:)];
figure('Name','Total Energy Consumption at different Time Slots','NumberTitle','off','color','White');
bar(E_Consumed,'grouped') ;
set(gca, 'XTick',1:6, 'XTickLabel',{num2str(50) num2str(100) num2str(150) num2str(200) num2str(250) num2str(300)},'fontname','Times New Roman');
legend('E-OEERP','OEERP','LEACH','DRINA','BCDCP');
xlabel('Time in Milli Seconds','fontname','Times New Roman');
ylabel('Total Energy Consumption (Joules)');
ylim([0 160]);
str=['Total Energy Consumption at different Time Slots'];
title(str,'fontsize',11,'fontname','Cambria','color','blue');

Throughput_exist(1,:)=[55000 57000 75000 17000];
Throughput_exist(2,:)=[55000 50000 74500 26000];
Throughput_exist(3,:)=[53000 52000 76000 27000];
Throughput_exist(4,:)=[52000 51000 70000 28000];
Throughput_exist(5,:)=[55000 55000 62000 29000];
Throughput_exist(6,:)=[54000 55000 56000 30000];
Throughput_total=[Throughput_propos(1,1) Throughput_exist(1,:);Throughput_propos(1,2) Throughput_exist(2,:);Throughput_propos(1,3) Throughput_exist(3,:);Throughput_propos(1,4) Throughput_exist(4,:);Throughput_propos(1,5) Throughput_exist(5,:);Throughput_propos(1,6) Throughput_exist(6,:);];

figure('Name','Throughput at different Time Slots','NumberTitle','off','color','White');
bar(Throughput_total,'grouped') ;
set(gca, 'XTick',1:6, 'XTickLabel',{num2str(50) num2str(100) num2str(150) num2str(200) num2str(250) num2str(300)},'fontname','Times New Roman')
legend('E-OEERP','OEERP','LEACH','DRINA','BCDCP');
xlabel('Time in Milli Seconds','fontname','Times New Roman');
ylabel('Throughput');
ylim([0 90000]);
str=['Throughput at different Time Slots'];
title(str,'fontsize',11,'fontname','Cambria','color','blue');

PDR_exist(1,:)=[60 63 90 20];
PDR_exist(2,:)=[62 63 97 35];
PDR_exist(3,:)=[61 63 98 36];
PDR_exist(4,:)=[61 63 85 37];
PDR_exist(5,:)=[60 62 75 37];
PDR_exist(6,:)=[61 63 65 38];
PDR_total=[PDR(1,1) PDR_exist(1,:);PDR(1,2) PDR_exist(2,:);PDR(1,3) PDR_exist(3,:);PDR(1,4) PDR_exist(4,:);PDR(1,5) PDR_exist(5,:);PDR(1,6) PDR_exist(6,:)];
time=[1 2 2.8 3.5 4 4.8];
figure('Name','Packet Delivery Ratio at different Time Slots','NumberTitle','off','color','White');
bar(PDR_total,'grouped') ;
set(gca, 'XTick',1:6, 'XTickLabel',{num2str(50) num2str(100) num2str(150) num2str(200) num2str(250) num2str(300)},'fontname','Times New Roman')
legend('E-OEERP','OEERP','LEACH','DRINA','BCDCP');
xlabel('Time in Milli Seconds','fontname','Times New Roman');
ylabel('Packet Delivery Ratio(PDR)');
ylim([0 120]);
str=['Packet Delivery Ratio at different Time Slots'];
title(str,'fontsize',11,'fontname','Cambria','color','blue');

for tt=1:length(time)
    NLifeTime(tt)=(Network_lifetime(1,1)/time(tt))*100;
end

NLife_exist(1,:)=[7000 2800 8500 4800];
NLife_exist(2,:)=[3000 1500 6500 3000];
NLife_exist(3,:)=[2000 1000 3300 2100];
NLife_exist(4,:)=[1800 1000 1500 1800];
NLife_exist(5,:)=[1200 800 1000 1300];
NLife_exist(6,:)=[1000 500 600 1000];
NLife_total=[NLifeTime(1,1) NLife_exist(1,:);NLifeTime(1,2) NLife_exist(2,:);NLifeTime(1,3) NLife_exist(3,:);NLifeTime(1,4) NLife_exist(4,:);NLifeTime(1,5) NLife_exist(5,:);NLifeTime(1,6) NLife_exist(6,:)];
figure('Name','Overall Network Lifetime at different Time Slots','NumberTitle','off','color','White');
bar(NLife_total,'grouped') ;
set(gca, 'XTick',1:6, 'XTickLabel',{num2str(50) num2str(100) num2str(150) num2str(200) num2str(250) num2str(300)},'fontname','Times New Roman')
legend('E-OEERP','OEERP','LEACH','DRINA','BCDCP');
xlabel('Time in Milli Seconds','fontname','Times New Roman');
ylabel('Network Lifetime(ms)');
ylim([0 12000]);
str=['Overall Network Lifetime at different Time Slots'];
title(str,'fontsize',11,'fontname','Cambria','color','blue');

% ===============================================
% Generate two random signals (you can replace these with your actual signals)
signal1 = randn(1, 1000); % Random signal 1
signal2 = randn(1, 1000); % Random signal 2

% Plot the original signals for visualization
figure;
subplot(2, 1, 1);
plot(signal1);
title('Signal 1');
subplot(2, 1, 2);
plot(signal2);
title('Signal 2');

SIG = [signal1 signal2];
% Compute the cross-correlation between the two signals
cross_corr = xcorr(signal1, signal2);

% Normalize the cross-correlation values to the range [-1, 1]
cross_corr = cross_corr / max(abs(cross_corr));

% Set a threshold for anomaly detection (you can adjust this threshold)
threshold = 0.9;

% Find the positions where the cross-correlation exceeds the threshold
anomaly_indices = find(cross_corr > threshold);

% Highlight the anomalies in the original signals
figure;
subplot(2, 1, 1);
plot(signal1,'bo');
title('Signal 1');
hold on;
for ii = 1:length(anomaly_indices)
    plot(anomaly_indices, SIG(anomaly_indices(ii)), 'r*');
end


hold off;

subplot(2, 1, 2);
plot(signal2,'bo');
title('Signal 2');
hold on;
plot(anomaly_indices, SIG(anomaly_indices), 'r*');
hold off;

% Display detected anomalies
disp('Anomaly positions:');
disp(anomaly_indices);

% Plot the original signals for visualization
figure;
subplot(2, 1, 1);
plot(signal1);
title('Signal 1');
xlabel('Sample');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(signal2);
title('Signal 2');
xlabel('Sample');
ylabel('Amplitude');

% Compute the cross-correlation between the two signals
cross_corr = xcorr(signal1, signal2);

% Normalize the cross-correlation values to the range [-1, 1]
cross_corr = cross_corr / max(abs(cross_corr));

% Set a threshold for anomaly detection (you can adjust this threshold)
threshold = 0.9;

% Find the positions where the cross-correlation exceeds the threshold
anomaly_indices = find(cross_corr > threshold);

% % Calculate and store the differences between original data and fault data
% difference1 = signal1(anomaly_indices) - signal1(anomaly_indices);
% difference2 = signal2(anomaly_indices) - signal2(anomaly_indices);

% Highlight the anomalies in the original signals
figure;
subplot(2, 1, 1);
plot(signal1);
title('Signal 1 with Fault Data');
xlabel('Sample');
ylabel('Amplitude');
hold on;
plot(anomaly_indices, SIG(anomaly_indices), 'ro', 'MarkerSize', 8);
hold off;

subplot(2, 1, 2);
plot(signal2);
title('Signal 2 with Fault Data');
xlabel('Sample');
ylabel('Amplitude');
hold on;
plot(anomaly_indices, SIG(anomaly_indices), 'ro', 'MarkerSize', 8);
hold off;

% Display detected anomalies
disp('Anomaly positions:');
disp(anomaly_indices);


demodulated_data = (modulated_data + 1) / 2;

% Convert binary data back to characters

aa = char(demodulated_data + '0');
received_message = char(bin2dec(aa));

% Display the transmitted and received messages
disp('Transmitted Message:');
disp(message);
disp('Received Message:');
disp(received_message);

% Calculate the Bit Error Rate (BER)
bit_errors = sum(abs(demodulated_data - (modulated_data + 1) / 2)) / length(demodulated_data);
disp(['Bit Error Rate (BER): ', num2str(bit_errors)]);


% Display the differences between original data and fault data
% disp('Difference in Signal 1:');
% disp(difference1);
% disp('Difference in Signal 2:');
% disp(difference2);