%% Main script to operate the thorlab spectrophotometer


addpath('./HardwareControlModule/');
addpath('./ADEModule/');
addpath('./IntegralControlModule/');
addpath('./VisualizationModule/');


% Initialize the spectrophotometer
TS = ThorSpectrophotometer();
TS.inttime = 0.00012;

h = msgbox('End Program');

while(ishandle(h))
tic
[wldata,specdata] = TS.getSpectrum;
plot(wldata.value,specdata.value);
ylim([0 1]);
toc
drawnow
pause(0.2);
end 

% tic
% time  = zeros(1,10);
% spec430 = zeros(1,10);
% spec630 = zeros(1,10);
% 
% for i = 1:10
%     binData = TS.getBinnedData([430, 630], 2.5)
%     time(i) = toc;
%     spec430(i) = binData(1);
%     spec630(i) = binData(2);
%     pause(1);
% end
% plot(time, spec630);
% hold on
% plot(time, spec430);

TS.delete();
