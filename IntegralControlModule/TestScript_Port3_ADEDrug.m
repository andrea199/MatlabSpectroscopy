%% Script to test the integral module
% In the future: Part of this script will be integrated into 
% the visualization class.
addpath('../ADEModule/');
addpath('../VisualizationModule/');
addpath('../HardwareControlModule');

pumpsFlag = 1;

%Initialize pumps
if pumpsFlag
KDc =KDPump('COM4');
KDd =KDPump('COM7');
end

%User inputs: 
tmax = 800;
dd =   3*(1e-6/3600); %m^3/s
qcss = 10*((1e-6)/3600); 
qcsm = 30*(1e-6/3600); 
cd  = 1;
numPumps = 2;
numDrugs = 1;

%Obtain Radius as a function of lenth
inrad = 2.1/1000; %m
Ucin = (qcsm+dd)/(pi*inrad^2);
nGrid= 0.5;
r1 = ones(1,length(0:nGrid:100-nGrid))*inrad; %More accurate to use 1.3 instead of 1.5; (changed 80 to 108)
r2 = ones(1,length(100:nGrid:204-nGrid))*0.8/1000;        %Distance measured is 104
r3 = ones(1,length(204:nGrid:390))*0.6/1000;           %Distance measured is 190
Rtot = [r1 r2 r3];
Rtots = smooth(Rtot,1);               

%ADE Solver
ADE = ADEModule();
ADE.Ngrid = length(Rtots);
ADE.RoLScale = Rtots;
ADE.UScale = Ucin*(inrad^2./Rtots.^2);
ADE.DeltaT = 1;
ADE.DomainLength = 0.390;
D = 1e-9;
ADE.PecletNum =  ADE.UScale.*ADE.RoLScale / D;
ADE.Dcoeff =  D*(1 + (ADE.PecletNum.^2)/48);
Nt = 200;
ADE.initializeCrankNicolson();


%Initialize the ADEModuleDrug
ADEDrug = ADEModule();
ADEDrug.Ngrid =  20;
ADEDrug.RoLScale = ones(ADEDrug.Ngrid,1)*inrad;
Udcin = dd/(pi*inrad^2);
ADEDrug.UScale = ones(ADEDrug.Ngrid,1)*Udcin;
ADEDrug.DeltaT = 1;
ADEDrug.DomainLength = 0.015;
D = 1e-9;
ADEDrug.PecletNum =  ADEDrug.UScale.*ADEDrug.RoLScale / D;
ADEDrug.Dcoeff =  D*(1 + (ADEDrug.PecletNum.^2)/48);
ADEDrug.initializeCrankNicolson();




%Initalize the Integral Control Module
ICM = IntegralControlModule(ADE);
ICM.initializePumps(numPumps, [qcsm,0.1*qcsm]);
ICM.initalizeDrugs(numDrugs, [qcss, dd/cd], cd);

%Initalize the Visualization Module
VM = VisualizationModule(ICM);

%Graphs:
close all
n1 = 1;
filename = 'pump.ai_ParkerAlgorithmRealTime.gif'; %name of gif

%Log time history
CarrierFlowRate = zeros (1,length(1:tmax/ADE.DeltaT));
DrugFlowRate = CarrierFlowRate;

i = 1;
sysTime = tic;
tic
pause(1);

for ind=1:tmax/ADE.DeltaT
    t=(ind-1)*ADE.DeltaT;
    %Drug
    ADEDrug.stepCrankNicolson();
    qdrug = ADEDrug.Solver.q;
    
    %Carrier/Mixed Stream
    ADE.Solver.q(1)= qdrug(end);
    ADE.stepCrankNicolson();
    q = ADE.Solver.q;
    ADE.Solver.qcorr = ADE.Solver.q.*(inrad^2./Rtots.^2);
    ICM.nextStepPrediction(ADE);
    
    % Plotting numerical solution
    VM.visualizeDrugDelivery(ADE,ICM,t); 
    
    % Store values
    CarrierFlowRate(i) = ICM.Pumps.PumpFlowRate(1);
    DrugFlowRate(i)   = ICM.Pumps.PumpFlowRate(2);
    
    
    %Update matrices for drug
    Udin = (DrugFlowRate(i))/(pi*inrad^2);
    ADEDrug.UScale = ones(ADEDrug.Ngrid,1)*Udcin;
    ADEDrug.PecletNum =  ADEDrug.UScale.*ADEDrug.RoLScale / D;
    ADEDrug.Dcoeff =  D*(1 + (ADEDrug.PecletNum.^2)/48);
    ADEDrug.initializeCrankNicolson();
    ADEDrug.Solver.q = qdrug;
    
    
    %Update matrices for carrier
    Ucin = (CarrierFlowRate(i)+DrugFlowRate(i))/(pi*inrad^2);
    ADE.UScale = Ucin*(inrad^2./Rtots.^2);
    ADE.PecletNum =  ADE.UScale.*ADE.RoLScale / D;
    ADE.Dcoeff =  D*(1 + (ADE.PecletNum.^2)/48);
    ADE.initializeCrankNicolson();
    ADE.Solver.q = q;
    
     if toc(sysTime)-t <1
         pause(1);
     end
    
    % Run pumps:  
    if pumpsFlag
    if toc>1
        diff=toc(sysTime)-t-1
        try
         KDc.infuseRate = (ICM.Pumps.PumpFlowRate(1)/(1e-6/3600)); % Flow rate set in ml/hr
        catch
        end
        try
        KDd.infuseRate = (ICM.Pumps.PumpFlowRate(2)/(1e-6/3600)); % Flow rate set in ml/hr
        catch  
        end
    tic
    end
    end
       
    
    i = i+1;
   % create file in gif
    frame = getframe(1);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    if n1 == 1
 		imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.05);
        n1=0;
 	else
		imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.05);
    end
end
