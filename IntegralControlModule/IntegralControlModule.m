classdef IntegralControlModule < handle
    %Class with routines to control syringe pumps through an integral
    %control approach
    
    properties
        UMax = []; % Maximum flow rate in the system -> updated from GUI (currently cast as velocity)
        UMaxfactor = 10; %Ratio of max flow rate to steady state flow rate
        Pumps; % Structure to store pump data
        Drugs; % Structure to store drug data
        DrugDeliveryHistoryTarget = 0; % Drug delivery historical target
        DrugDeliveryHistoryCurrent = 0;% Drug delivery cummulative value
        conFactor; %Convert ml/hr to m^3/s
        %Multi Drug Control properties 
        Wt; %Weighting factor
        ConcRange = 0.05; %Conc fluctuation range
        cd
        Ci;
        Qci;
        Qcmax;
        Cimax;
        ConcEq; % Concentration at equilibrium
        fun;
        A;
        b;
        Aeq;
        beq;
        lb;
        ub;
        options;
    end
    
   methods % Constructor
        function obj = IntegralControlModule()
            % Removed the need to input the ADE object
        end     
   end 
     
   methods % Core functions:
        function nextStepPrediction(obj,ADEModuleObj)
            % Predicts the pump flow rates for the next time step
            
            % Catch error:
            if(~isfield(obj.Pumps,'NumPumps'))
                error('Run initializePumps function in the IntegralControlModule');
            end
            
            if(~isfield(obj.Drugs,'NumDrugs'))
                error('Run initalizeDrugs function in the IntegralControlModule');
            end
            
            % Evolve using Parkers Algorithm
            %ParkersAlgorithm(obj,ADEModuleObj)
            %IntegralControlAlgorithm(obj,ADEModuleObj);   
            multiDrugControlAlgorithm(obj,ADEModuleObj);   
            
        end
        
        function ParkersAlgorithm(obj,ADEModuleObj)
            % Implemention of Parkers algorithm:
            obj.Pumps.TotalPumpFlowRate  = abs(sum(obj.Drugs.DrugFlowRate)/(ADEModuleObj.Solver.qcorr(end)));                            
            if obj.Pumps.TotalPumpFlowRate > sum(obj.Pumps.MaxPumpFlowRate)
                obj.Pumps.PumpFlowRate(1) = obj.Pumps.MaxPumpFlowRate(1);
                obj.Pumps.PumpFlowRate(2) = obj.Pumps.PumpFlowRate(1)* (obj.Drugs.DrugFlowRate(2)/obj.Drugs.DrugFlowRate(1));
            else
                obj.Pumps.PumpFlowRate(1) = obj.Pumps.TotalPumpFlowRate/(1 +  (obj.Drugs.DrugFlowRate(2)/obj.Drugs.DrugFlowRate(1)));
                obj.Pumps.PumpFlowRate(2) = obj.Pumps.TotalPumpFlowRate - obj.Pumps.PumpFlowRate(1);
            end
        end
        
        function IntegralControlAlgorithm(obj,ADEModuleObj)
            % Implementation of Time History (Integral) Control Algorithm
            
            %Expected cumulative value:
            obj.DrugDeliveryHistoryTarget = obj.DrugDeliveryHistoryTarget + obj.Drugs.DrugFlowRate(2)*obj.Drugs.DrugStockConc*ADEModuleObj.DeltaT;
            obj.DrugDeliveryHistoryCurrent = obj.DrugDeliveryHistoryCurrent + ADEModuleObj.Solver.q(end)*sum(obj.Pumps.PumpFlowRate)*(obj.Drugs.DrugFlowRate(2)/sum(obj.Drugs.DrugFlowRate))*obj.Drugs.DrugStockConc*ADEModuleObj.DeltaT;

            %Current value:
            obj.Pumps.TotalPumpFlowRate  = abs(sum(obj.Drugs.DrugFlowRate)/(ADEModuleObj.Solver.q(end)));     
            if obj.Pumps.TotalPumpFlowRate > sum(obj.Pumps.MaxPumpFlowRate)||  obj.DrugDeliveryHistoryTarget > obj.DrugDeliveryHistoryCurrent
                obj.Pumps.PumpFlowRate(1) = obj.Pumps.MaxPumpFlowRate(1);
                obj.Pumps.PumpFlowRate(2) = obj.Pumps.PumpFlowRate(1)* (obj.Drugs.DrugFlowRate(2)/obj.Drugs.DrugFlowRate(1));
            else
                obj.Pumps.PumpFlowRate(1) = obj.Pumps.TotalPumpFlowRate/(1 +  (obj.Drugs.DrugFlowRate(2)/obj.Drugs.DrugFlowRate(1)));
                obj.Pumps.PumpFlowRate(2) = obj.Pumps.TotalPumpFlowRate - obj.Pumps.PumpFlowRate(1);
            end
        end
        
        
        function multiDrugControlAlgorithm(obj, ADEModuleObj)
            % Nonlinear optimization based control algorithm (dev)
            
            %Current state
            obj.cd = ADEModuleObj.Solver.qcorr(end,:);
            %obj.Wt = 1;
            
            %Update constraint matrix
            obj.A(end-obj.Drugs.NumDrugs+1:end, 2:end)=eye(obj.Drugs.NumDrugs).*obj.cd'; 
   
            % Initial guess
            q0 = obj.Pumps.PumpFlowRate'/obj.conFactor; 
            % Solver
            q = fmincon(obj.fun,q0,obj.A,obj.b,obj.Aeq,obj.beq,obj.lb,obj.ub,[], obj.options);    
            
            %Assign flow rates 
            obj.Pumps.PumpFlowRate(1) = q(1)*obj.conFactor;
            obj.Pumps.PumpFlowRate(2:end) = q(2:end)*obj.conFactor;
        end
        
   end
   
   methods %Supporting functions
       function initializeMultiDrugControlAlgorithm(obj)  
           %Initialize the multi drug algorithm including
           % (a) Objective function
           % (b) Assemble constraint matrices
           % v.2 Multi drug
           
            %Initial values
            obj.conFactor = (1e-6/3600); % Convert flowrate to m^3/s
            obj.Qci = obj.Drugs.DrugFlowRate(1)/obj.conFactor; 
            obj.Ci = obj.Drugs.DrugFlowRate(2:end)/obj.conFactor; 
            obj.Qcmax = obj.Pumps.MaxPumpFlowRate(1)/obj.conFactor;
            obj.Cimax = obj.Pumps.MaxPumpFlowRate(2:end)/obj.conFactor;
            obj.ConcEq =  obj.Ci./(obj.Ci+obj.Qci);
                        
            % Constraints            
            obj.Pumps.ModMaxCarrierPumpFlowRate = obj.Qci + (1-obj.Wt)*(obj.Qcmax - obj.Qci);
%              obj.A = [eye(obj.Pumps.NumPumps); -eye(obj.Pumps.NumPumps); [ zeros(obj.Drugs.NumDrugs,1) eye(obj.Drugs.NumDrugs)]]; % Generalized for multiple drugs
%              obj.b =[[obj.Pumps.ModMaxCarrierPumpFlowRate; obj.Pumps.MaxPumpFlowRate(2:end)'/obj.conFactor]; [0;-0.5*obj.Drugs.DrugFlowRate(2:end)'/obj.conFactor];1.1*obj.Ci']; %Uses modified max flowrates.
%              obj.Aeq = [-ones(obj.Drugs.NumDrugs,1).*(obj.ConcEq./(1-obj.ConcEq))' eye(obj.Drugs.NumDrugs)];
%              obj.beq =[zeros(obj.Drugs.NumDrugs,1)];
             obj.A = [eye(obj.Pumps.NumPumps); -eye(obj.Pumps.NumPumps); [-ones(obj.Drugs.NumDrugs,1).*((1+obj.ConcRange)*obj.ConcEq./(1-(1+obj.ConcRange)*obj.ConcEq))' eye(obj.Drugs.NumDrugs)]; [ones(obj.Drugs.NumDrugs,1).*((1-obj.ConcRange)*obj.ConcEq./(1-(1-obj.ConcRange)*obj.ConcEq))' -eye(obj.Drugs.NumDrugs)]; [ zeros(obj.Drugs.NumDrugs,1) eye(obj.Drugs.NumDrugs)]]; % Generalized for multiple drugs
             obj.b =[[obj.Pumps.ModMaxCarrierPumpFlowRate; obj.Pumps.MaxPumpFlowRate(2:end)'/obj.conFactor]; [0;-0.5*obj.Drugs.DrugFlowRate(2:end)'/obj.conFactor];zeros(obj.Drugs.NumDrugs,1);zeros(obj.Drugs.NumDrugs,1);1.5*obj.Ci'];
             obj.Aeq = [];
             obj.beq = [];
             obj.lb = [];
             obj.ub = [];
            
            obj.options = optimoptions('fmincon','Display','iter','Algorithm','sqp'); %Sqp alogrith appears to be more stable

       end
      
       function initializeObjectiveFunction(obj)
           %Initializes the objective function based on the number of drugs in the system.
           %General form of the control algorithm:  obj.fun = @(q) exp(abs(obj.Ci - (q(2))*cd)) + (obj.Wt*(obj.Qci - q(1))^2 + (1-obj.Wt)*(1-cd)*(obj.Qcmax - q(1)).^2)/(obj.Qcmax - obj.Qci); 
           switch obj.Drugs.NumDrugs
               case 1
%              obj.fun = @(q) exp(abs(obj.Ci - (q(2))*obj.cd)) + (obj.Wt*(obj.Qci - q(1))^2 + (1-obj.Wt)*(1-obj.cd)*(obj.Qcmax - q(1)).^2)/(obj.Qcmax - obj.Qci); 
%                obj.fun = @(q) exp(abs(obj.Ci - (q(2))*obj.cd)) + ((1+obj.cd)^2*(obj.Wt)*(obj.Qci - q(1))^2 + (1-obj.Wt)*(1-obj.cd)^2*(obj.Qcmax - q(1)).^2 )/(obj.Qcmax - obj.Qci);            
               obj.fun = @(q) exp(abs(obj.Ci - (q(2))*obj.cd)) + ((1+obj.cd)^2*(obj.Wt)*(obj.Qci - q(1))^2 + (1-obj.Wt)*(max([(1-obj.cd),0])^2*(obj.Qcmax - q(1)).^2 + ...
               max([(obj.cd-1),0])^2*(0-q(1)).^2))/(obj.Qcmax - obj.Qci) + (obj.cd)*(1+obj.Wt)*(abs(obj.Pumps.PumpFlowRate(2)/obj.conFactor -q(2)))/(obj.Drugs.DrugFlowRate(2)/obj.conFactor);  %0.01 for simple geometries    % Might be required to prevent overshoots when the concentration constraint is relaxed.  % Damping term added to prevent oscillations     

               case 2 
%                    obj.fun = @(q) exp(abs(obj.Ci(1) - (q(2))*obj.cd(1))) + exp(abs(obj.Ci(2) - (q(3))*obj.cd(2))) + ...
%                    (obj.Wt*(obj.Qci - q(1))^2 + (1-obj.Wt)*(1-obj.cd(1))*(obj.Qcmax - q(1)).^2)/(obj.Qcmax - obj.Qci) + ...
%                    (obj.Wt*(obj.Qci - q(1))^2 + (1-obj.Wt)*(1-obj.cd(2))*(obj.Qcmax - q(1)).^2)/(obj.Qcmax - obj.Qci);   % Base case
%             
%                   obj.fun = @(q) exp(abs(obj.Ci(1) - (q(2))*obj.cd(1))) + exp(abs(obj.Ci(2) - (q(3))*obj.cd(2))) + ...
%                   (obj.Wt*(1+obj.cd(1))^2*(obj.Qci - q(1))^2 + (1-obj.Wt)*(1-obj.cd(1))^2*(obj.Qcmax - q(1)).^2)/(obj.Qcmax - obj.Qci) + ...
%                   (obj.Wt*(1+obj.cd(2))^2*(obj.Qci - q(1))^2 + (1-obj.Wt)*(1-obj.cd(2))^2*(obj.Qcmax - q(1)).^2)/(obj.Qcmax - obj.Qci);   % (1+cd)^2

                    obj.fun = @(q) exp(abs(obj.Ci(1) - (q(2))*obj.cd(1))) + exp(abs(obj.Ci(2) - (q(3))*obj.cd(2))) + ...
                    (obj.Wt*(1+obj.cd(1))^2*(obj.Qci - q(1))^2 + (1-obj.Wt)*max([(1-obj.cd(1)),0])^2*(obj.Qcmax - q(1)).^2  + (1-obj.Wt)*max([(obj.cd(1)-1),0])^2*(0 - q(1)).^2)/(obj.Qcmax - obj.Qci) + ...
                    (obj.Wt*(1+obj.cd(2))^2*(obj.Qci - q(1))^2 + (1-obj.Wt)*max([(1-obj.cd(2)),0])^2*(obj.Qcmax - q(1)).^2 + (1-obj.Wt)*max([(obj.cd(2)-1),0])^2*(0 - q(1)).^2)/(obj.Qcmax - obj.Qci);   % (1+cd)^2 and accounting for overshooting 


           end
       end   
       
      function initializePumps(obj,NumPumps,MaxFlowRate)
          % Intializes the number of active pumps in the system
          % Inputs: NumPumps - Number of active pumps
          %         MaxFlowRate - Maximum allowed flow rate for the systems
          %         (The input format is an array (1,NumPumps))
                   
          obj.Pumps.NumPumps = NumPumps;
          obj.Pumps.PumpFlowRate = MaxFlowRate; %Initial flow rates set at their max vals
          obj.Pumps.MaxPumpFlowRate = MaxFlowRate;

      end
      % function addPump(obj)
      
      function initalizeDrugs(obj, NumDrugs, DrugFlowRate,DrugStockConc)
         % Initializes the number of active drugs in the system
         % Inputs: NumDrugs - Number of active pumps (Drug 1 is carrier)
         %         DrugFlowRate - Intended mass flow rate for the drugs at
         %         stead state. (The input format is a array (1,NumDrugs) )
         %         DrugStockConc - The stock concentration of the drug
         
         obj.Drugs.NumDrugs = NumDrugs;
         obj.Drugs.DrugFlowRate = DrugFlowRate; 
         obj.Drugs.DrugStockConc = DrugStockConc;       
      end 
      
      function resetPumps(obj)
        % Resets pump flow rates to initial values
        obj.Pumps.PumpFlowRate(1) = obj.Pumps.ModMaxCarrierPumpFlowRate*obj.conFactor;
        %obj.Pumps.PumpFlowRate(2:end) = obj.Pumps.ModMaxCarrierPumpFlowRate * obj.ConcEq/(1 - obj.ConcEq)*obj.conFactor; %Not compatible with the objective function   
      end
      
  end
    
end