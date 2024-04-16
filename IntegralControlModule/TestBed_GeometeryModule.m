%% Test bed to test the performance of the control law

%Link Libraries

addpath('../ADEModule/');
addpath('../VisualizationModule/');
addpath('../HardwareControlModule/');
addpath('../GeometryModule/');

%User parameters:
tmax = 1000; 
ports = [3]; %Remove/Add drug 2
numDrugs = length(ports);
numPumps = numDrugs+1;

qDrugSteady = [3]*(1e-6/3600); qDrugMax = [10]*(1e-6/3600);  %max 10
qDrug2Steady = [5]*(1e-6/3600); qDrug2Max = [15]*(1e-6/3600); %max 15, min 5
qSalineSteady = [10]*(1e-6/3600); qSalineMax = [30]*(1e-6/3600);%max 30

%Drug Parameters
drugDiff = 1e-9;

%Geometry 
CG = CatherGeometry('BBN', '456003');
CG.parseGeometry(ports, 0.5); %first argument is the port number (port 1 is closest to patient), second argument is nGrid

%Initializ the ICM
ICM = IntegralControlModule(); 
ICM.initializePumps(numPumps, [qSalineMax,qDrugMax]);
ICM.initalizeDrugs(numDrugs, [qSalineSteady, qDrugSteady], 1); 
ICM.Wt = 0.05;
ICM.ConcRange = 0.25;
ICM.initializeMultiDrugControlAlgorithm();
ICM.initializeObjectiveFunction();

%Initalize the Visualization Module
VM = VisualizationModule(ICM);

%Initialize the MultiDrugSolver
MDS = MultiDrugSolver(CG,ICM);
MDS.drugDiffusivity=drugDiff;

%Initialize var to store data
time = zeros(1, 1000);
psi  = time;
psi2 = time;
salinePump = time;
drugPump   = time;

%Evolve over time
 MDS.CalculateFirstTimeStep();
 
 tic
for ind=1:tmax/MDS.DeltaT
    t=(ind-1)*MDS.DeltaT;
    
    %Store data
    time(ind) = t;
    psi(ind) = MDS.ADE{1}.Solver.qcorr(end,1);
    %psi2(ind) = MDS.ADE{1}.Solver.qcorr(end,2);
    salinePump(ind) = ICM.Pumps.PumpFlowRate(1);
    drugPump (ind) = ICM.Pumps.PumpFlowRate(2);
    %Plotting numerical solution
    %VM.visualizeDrugDelivery(MDS.ADE{1},ICM,t); 
   
    %Next step prediction
    ICM.nextStepPrediction(MDS.ADE{1});
       
    %Calculate the next step
    MDS.CalculateNextTimeStep();
end
toc
% Optionally write to text file
FineSalineInletArea = 1.40992e-05;
FineDrugPort3Area = 1.44383e-05;
FineDrugPort1Area =1.40778e-05;

formatSpec = '(%03d (%7.4e 0 0))\n';
fileID = fopen('InletVelocity.txt','w');
fprintf(fileID,formatSpec, [time;salinePump./FineSalineInletArea]);
fclose(fileID);
fileID = fopen('Port3Velocity.txt','w');
fprintf(fileID,formatSpec, [time;-drugPump./FineDrugPort3Area]);
fclose(fileID);


%To do:
% Run the simulation in the cluster
% Scale of the flow rates.
% (i)Check if the intialization with max pump flow rates has any effect on stability?
% (ii) Calling  VM.visualizeDrugDelivery(MDS.ADE{1},ICM,t);  before
% ICM.nextStep prediction had a very different plotted result. 
% (iii)Port 3 numerical calculations were slighlty off (the initial filling was not appropriate) 