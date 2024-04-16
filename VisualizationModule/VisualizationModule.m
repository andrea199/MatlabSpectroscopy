classdef VisualizationModule < handle
    %Class with routines to visualize drug delivery
    
    properties
        PlotHandle = []; %Handle to the plot window
        PlotAxis; %Structure with handles to the axess
        DrugConcPlotHandle = [0 0 0 0]; % Handle to drug concentration
        DrugDeliveredHistory = []; % Handle to store the drug delivery history 
        TMax = 10; % Default Tmax for visualization plots
        TInjection = []; %Handle containing the injection history of the drugs
        ColorPanel = {[0.33 0.33 0.33],[192,192,255]/255, [255,208,160]/255, [128,255,128]/255}; % Color panel for graphs
        ColorPanelIntended = {'k','b','r','g'}; % Color panel for graph limits
    end
    
   methods % Constructor
        function obj = VisualizationModule(ICModuleObj)            
            obj.DrugDeliveredHistory = zeros(1,ICModuleObj.Drugs.NumDrugs);
            obj.TInjection = zeros(1,ICModuleObj.Drugs.NumDrugs);
        end     
   end 
     
   methods % Core functions:
        function visualizeDrugDelivery(obj,ADEModuleObj, ICModuleObj,time)
            % Plots the solution, drug delivery, and pump inputs
          
            % Check if plot window is available
            if (~sum(ishandle(obj.PlotHandle)))
                obj.initializePlotWindows(ADEModuleObj,ICModuleObj);
            else
                plot(obj.PlotAxis.l2,time,ICModuleObj.Pumps.PumpFlowRate(1),'.','Color',obj.ColorPanel{1},'Linewidth',1.4,'MarkerSize',20);   %Saline pump flow
                for i = 1:ICModuleObj.Drugs.NumDrugs
                     try
                     delete(obj.DrugConcPlotHandle(i))
                     catch
                     end
                     obj.DrugConcPlotHandle(i) = plot(obj.PlotAxis.l1,ADEModuleObj.Solver.x(1:ADEModuleObj.RealDomainIndex),ADEModuleObj.Solver.qcorr(:,i),'Color',obj.ColorPanel{i+1}, 'Linewidth',2.8,'MarkerSize',20);
                    
                     plot(obj.PlotAxis.l2,time,ICModuleObj.Pumps.PumpFlowRate(i+1),'.','Color',obj.ColorPanel{i+1},'Linewidth',1.4,'MarkerSize',20);  %Drug pump flow
                     obj.DrugDeliveredHistory(i) = obj.DrugDeliveredHistory(i) + ADEModuleObj.Solver.qcorr(end,i)*ICModuleObj.Pumps.PumpFlowRate(i+1)*ICModuleObj.Drugs.DrugStockConc*ADEModuleObj.DeltaT;
                     plot(obj.PlotAxis.l3,time, obj.DrugDeliveredHistory(i),'.','Color',obj.ColorPanel{i+1},'Linewidth',1.4,'MarkerSize',20);  
                     plot(obj.PlotAxis.l4,time, ADEModuleObj.Solver.qcorr(end,i)*ICModuleObj.Pumps.PumpFlowRate(i+1)*ICModuleObj.Drugs.DrugStockConc,'.','Color',obj.ColorPanel{i+1},'Linewidth',1.4,'MarkerSize',20); 
                end
                
                drawnow
                
                %Check plot range:
                if obj.TMax < time
                    obj.TMax = 10* obj.TMax;
                    obj.updatePlotAxisRange(ADEModuleObj,ICModuleObj);
                    obj.plotAnticipatedProfiles(ADEModuleObj,ICModuleObj);
                end
            end
            
             
        end
   end
   
   methods %Supporting functions
      function initializePlotWindows(obj,ADEModuleObj,ICModuleObj,plotAxisHandle)
          % Handle optional arguments
            if nargin>3
              obj.PlotHandle = plotAxisHandle;
              figure(plotAxisHandle);
            else  
              obj.PlotHandle = figure(1); 
              set(obj.PlotHandle, 'Position', [1 41 1536 748.8000],'color','w','NumberTitle', 'off','Name', 'Pump.ai');            
            end
          
          % Initalizes the plot windows
            Factor = 1.25;
            obj.PlotAxis.l1=axes('Position',[0.1300    0.1100    0.3750    0.3366]);
            set(obj.PlotAxis.l1,'FontName','Helvetica','Fontsize',13*Factor,'Linewidth',1.1);
            ylabel(obj.PlotAxis.l1,{'c/c_0'},'Fontsize',13*Factor);
            xlabel(obj.PlotAxis.l1,{'Length'},'Fontsize',13*Factor);
            axis(obj.PlotAxis.l1,[0,ADEModuleObj.DomainLength,0,2]);
            hold on

            obj.PlotAxis.l2=axes('Position',[0.1300    0.6100    0.3750    0.3366]);
            set(obj.PlotAxis.l2,'FontName','Helvetica','Fontsize',13*Factor,'Linewidth',1.1);
            ylabel(obj.PlotAxis.l2,{'Pump Flow Rate'},'Fontsize',13*Factor);
            xlabel(obj.PlotAxis.l2,{'Time'},'Fontsize',13*Factor);
            axis(obj.PlotAxis.l2,[0,obj.TMax,0,ICModuleObj.Pumps.MaxPumpFlowRate(1)]);
            hold on

            obj.PlotAxis.l3=axes('Position',[0.5800    0.1100    0.3750    0.3366]);
            set(obj.PlotAxis.l3,'FontName','Helvetica','Fontsize',13*Factor,'Linewidth',1.1);
            ylabel(obj.PlotAxis.l3,{'Cumulative Drug delivered'},'Fontsize',13*Factor);
            xlabel(obj.PlotAxis.l3,{'Time'},'Fontsize',13*Factor);
            axis(obj.PlotAxis.l3,[0,obj.TMax,0,ICModuleObj.Drugs.DrugFlowRate(2)*ICModuleObj.Drugs.DrugStockConc*obj.TMax]);
            hold on


            obj.PlotAxis.l4=axes('Position',[0.5800    0.6100    0.3750    0.3366]);
            set(obj.PlotAxis.l4,'FontName','Helvetica','Fontsize',13*Factor,'Linewidth',1.1);
            ylabel(obj.PlotAxis.l4,{'Drug Delivery Rate'},'Fontsize',13*Factor);
            xlabel(obj.PlotAxis.l4,{'Time'},'Fontsize',13*Factor);
            axis(obj.PlotAxis.l4,[0,obj.TMax,0,3*ICModuleObj.Drugs.DrugFlowRate(2)*ICModuleObj.Drugs.DrugStockConc]);
           
            hold on
            
            % Plot anticipated profiles
            obj.plotAnticipatedProfiles(ADEModuleObj,ICModuleObj);
          
      end   
      function plotAnticipatedProfiles(obj,ADEModuleObj,ICModuleObj)
          % This function plots the anticipated profiles 
          
          %Plot the anticipated pump profiles:
          plot(obj.PlotAxis.l2, [0,obj.TMax], [ICModuleObj.Drugs.DrugFlowRate(1), ICModuleObj.Drugs.DrugFlowRate(1)],'--','Color',obj.ColorPanelIntended{1}, 'Linewidth',1.4,'MarkerSize',20); 
          
          for  i= 1:ICModuleObj.Drugs.NumDrugs
          plot(obj.PlotAxis.l2, [0,obj.TMax], [ICModuleObj.Drugs.DrugFlowRate(i+1), ICModuleObj.Drugs.DrugFlowRate(i+1)],'--','Color',obj.ColorPanelIntended{i+1}, 'Linewidth',1.4,'MarkerSize',20);
          
          %Plot the anticipated cumulative profiles: 
          plot(obj.PlotAxis.l3, [0,obj.TMax], [0,ICModuleObj.Drugs.DrugFlowRate(i+1)*ICModuleObj.Drugs.DrugStockConc*obj.TMax],'--','Color',obj.ColorPanelIntended{i+1}, 'Linewidth',1.4,'MarkerSize',20);

          %Plot the anticipated instanteous profiles: 
          plot(obj.PlotAxis.l4, [0,obj.TInjection(i),obj.TMax], [0,ICModuleObj.Drugs.DrugFlowRate(i+1)*ICModuleObj.Drugs.DrugStockConc, ICModuleObj.Drugs.DrugFlowRate(i+1)*ICModuleObj.Drugs.DrugStockConc],'--','Color',obj.ColorPanelIntended{i+1},'Linewidth',1.4,'MarkerSize',20);       
          end
      end
   
    
   function updatePlotAxisRange(obj,ADEModuleObj,ICModuleObj)
         % This function updates the plot axis range
         
          axis(obj.PlotAxis.l2,[0,obj.TMax,0,ICModuleObj.Pumps.MaxPumpFlowRate(1)]);
          axis(obj.PlotAxis.l3,[0,obj.TMax,0,ICModuleObj.Drugs.DrugFlowRate(2)*ICModuleObj.Drugs.DrugStockConc*obj.TMax]);
          axis(obj.PlotAxis.l4,[0,obj.TMax,0,3*ICModuleObj.Drugs.DrugFlowRate(2)*ICModuleObj.Drugs.DrugStockConc]);
   end
   end
end