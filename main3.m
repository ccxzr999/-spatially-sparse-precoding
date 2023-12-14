clear
clear all
%snrValuesDB = -40:5:0; %in Figures
snrValuesDB = 0; %in Figures
snrValues = 10.^(snrValuesDB./10);
tryNumber = 100; 
numberDataStreams = 1:1:4;
%% Parameters
parameters_main3 = containers.Map('KeyType','char','ValueType','any');%������һ���µ�containers.Mapʵ������ָ���˼����ͣ�KeyType��Ϊ�ַ���char����ֵ���ͣ�ValueType��Ϊ�������ͣ�any����
parameters_main3("numberTransmitAntennas") = 64; % Number of transmit antennas parameters_main['key1'] = 'value1';  % ʹ��[]�����������ֵ��
parameters_main3("numberRecieveAntennas") = 16; % Number of receive antennas
parameters_main3("numberDataStreams") = 1; % Number of data streams
parameters_main3("numberRFchains") = 4; % Number of RF chains for precoding and combining
parameters_main3("numberCluster") = 8; % Number of clusters
parameters_main3("numberRayPercluster") = 10; % Number of rays per cluster
parameters_main3("angularSpread") = 7.5; % Angular spread of 7.5 degree
spectralEffOptimal = zeros(tryNumber,length(snrValues));
spectralEffHybrid = zeros(tryNumber,length(snrValues));
spectralEffBeam = zeros(tryNumber,length(snrValues));
for s = 1 :length(numberDataStreams)
    SNR = snrValues;
    parameters_main3("numberDataStreams") = numberDataStreams(s); % Number of data streams
    for i = 1:tryNumber
        channel = channel_generation(parameters_main3);
        tempObj = OptimalUnconstraint(SNR,channel);
        spectralEffOptimal(i,s) = tempObj.spectralEfficiency;
        tempObj = HybridSparsePrecoding(SNR,channel);
        spectralEffHybrid(i,s) = tempObj.spectralEfficiency;
    end
end

spectralEffOptimalSNR = mean(spectralEffOptimal,1); 
spectralEffHybridSNR = mean(spectralEffHybrid,1);
% figure(); 
hold on
%l1 = plot(numberDataStreams,spectralEffOptimalSNR,'-s','color',[0 0.5 0],'LineWidth',2.0,'MarkerSize',8.0);
%l2 = plot(numberDataStreams,spectralEffHybridSNR,'-o','Color',[0 0.45 0.74],'LineWidth',2.0,'MarkerSize',8.0);hold on;
color = rand(1,3);
l1 = plot(numberDataStreams,spectralEffOptimalSNR,'-s','Color',color,'LineWidth',2.0,'MarkerSize',8.0, 'DisplayName', sprintf("Optimal Uns. %dx%d, ",parameters_main3("numberTransmitAntennas"),parameters_main3("numberRecieveAntennas")) );
l2 = plot(numberDataStreams,spectralEffHybridSNR,'-o','Color',color,'LineWidth',2.0,'MarkerSize',8.0,'DisplayName', sprintf("Hybrid Comb. %dx%d, ",parameters_main3("numberTransmitAntennas"),parameters_main3("numberRecieveAntennas") ));

%legend([l1,l2],'Optimal unconstrained precoding N_s','Hybrid precoding and combining ','Location','northwest','FontSize',15);
xlabel("SNR (dB)",'FontSize', 20)
ylabel("Spectral Efficiency(bits/s/Hz)",'FontSize', 20)
%% v2
parameters_main3("numberTransmitAntennas") = 256; % Number of transmit antennas
parameters_main3("numberRecieveAntennas") = 16; % Number of receive antennas
spectralEffOptimal = zeros(tryNumber,length(numberDataStreams));
spectralEffHybrid = zeros(tryNumber,length(numberDataStreams));
for a = 1:length(numberDataStreams)
    parameters_main3("numberDataStreams") = numberDataStreams(a);
    for i = 1:tryNumber
        channel = channel_generation(parameters_main3);
        tempObj = OptimalUnconstraint(SNR,channel);
        spectralEffOptimal(i,a) = tempObj.spectralEfficiency;
        tempObj = HybridSparsePrecoding(SNR,channel);
        spectralEffHybrid(i,a) = tempObj.spectralEfficiency;
    end
end
% Averaging Tries
spectralEffOptimalSNR = mean(spectralEffOptimal,1); 
spectralEffHybridSNR = mean(spectralEffHybrid,1);
hold on
color = rand(1,3);
l3 = plot(numberDataStreams,spectralEffOptimalSNR,'-s','Color',color,'LineWidth',2.0,'MarkerSize',8.0, 'DisplayName', sprintf("Optimal Uns. %dx%d, ",parameters_main3("numberTransmitAntennas"),parameters_main3("numberRecieveAntennas")) );
l4 = plot(numberDataStreams,spectralEffHybridSNR,'-o','Color',color,'LineWidth',2.0,'MarkerSize',8.0,'DisplayName', sprintf("Hybrid Comb. %dx%d, ",parameters_main3("numberTransmitAntennas"),parameters_main3("numberRecieveAntennas") ));
%% v3
parameters_main3("numberTransmitAntennas") = 256; % Number of transmit antennas
parameters_main3("numberRecieveAntennas") = 64; % Number of receive antennas
spectralEffOptimal = zeros(tryNumber,length(numberDataStreams));
spectralEffHybrid = zeros(tryNumber,length(numberDataStreams));
for a = 1:length(numberDataStreams)
    parameters_main3("numberDataStreams") = numberDataStreams(a);
    for i = 1:tryNumber
        channel = channel_generation(parameters_main3);
        tempObj = OptimalUnconstraint(SNR,channel);
        spectralEffOptimal(i,a) = tempObj.spectralEfficiency;
        tempObj = HybridSparsePrecoding(SNR,channel);
        spectralEffHybrid(i,a) = tempObj.spectralEfficiency;
    end
end
% Averaging Tries
spectralEffOptimalSNR = mean(spectralEffOptimal,1); 
spectralEffHybridSNR = mean(spectralEffHybrid,1);
hold on
color = rand(1,3);
l5 = plot(numberDataStreams,spectralEffOptimalSNR,'-s','Color',color,'LineWidth',2.0,'MarkerSize',8.0, 'DisplayName', sprintf("Optimal Uns. %dx%d ",parameters_main3("numberTransmitAntennas"),parameters_main3("numberRecieveAntennas") ));
l6 = plot(numberDataStreams,spectralEffHybridSNR,'-o','Color',color,'LineWidth',2.0,'MarkerSize',8.0, 'DisplayName', sprintf("Hybrid Comb. %dx%d ",parameters_main3("numberTransmitAntennas"),parameters_main3("numberRecieveAntennas") ));
legend
xlabel("DataStreams (number)")
ylabel("Spectral Efficiency(bits/s/Hz)")