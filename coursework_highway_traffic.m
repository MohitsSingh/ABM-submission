% Agent-based modeling coursework
% Highway traffic simulation based on cellular automata
clc;
clear;

% This script contains the coded for both Phase Transition Model and
% Agent-based model.

load Full_Data1;% Full_Data1 contains all the real-time data for first 15mins
jamDens=0.26;% set up jam density
A=200;
B=500;

% setting the time scale
t_i=min(Full_Data1(:,3)/1000)+150; % find the minimum initial time in seconds
t_f=max(Full_Data1(:,3)/1000)-150; % find the maximum final time in seconds
                                   % and the first and last 150s are cut
% setting the Location scale
X_i=0+100;    % find the minimum initial Location in feet
X_f=1650-100; % find the maximum final Location in feet
             % and the first and last 100 feet are cut

% setting cell dimensions
t_length=4;  % giving that the cell dimension is 4second*70feet
X_length=70;

% setting cell array
t=t_i:t_length:t_f;dt=t(2)-t(1); % cells in horizontal dimension, dividing time in unit of time length
X=X_i:X_length:X_f;dX=X(2)-X(1); % cells in vertical axis, dividing Location in unit of Location length

%it is necessary to creat an zero matrix and then giving values for each element of this matrix
Measure_countArray=zeros(size(t,2)-1,size(X,2)-1); 

% creating an zeros matrix for number of vehicle in each cell
VehNum=zeros(size(Measure_countArray));
VehNum1=zeros(size(Measure_countArray));
VehNum2=zeros(size(Measure_countArray));

% creating zeros matrix for all necessary quantities
TotalVel=zeros(size(Measure_countArray));
TotalPertEst=zeros(size(Measure_countArray));
TotalDensEst=zeros(size(Measure_countArray));

%starting the data analysis 
[n_row,n_col]=size(Full_Data1);
index_begin=1;counter=0;
for i=1:n_row-1
    index_end=i;
    if Full_Data1(i,1)==Full_Data1(i+1,1)
        index_end=i;
    else
        SamVehicle=Full_Data1(index_begin:index_end,:);% collect all the rows for same vehicle
        if SamVehicle(1,7)>1&&SamVehicle(1,7)<6
            Time=SamVehicle(2:end-1,3)/1000;
            Location=SamVehicle(2:end-1,4);
            counter=counter+1;
            Velocity=SamVehicle(2:end-1,5);
            Acceleration=SamVehicle(2:end-1,6);
            for j=1:size(Time,1)
                if Time(j)>t(1)&&Time(j)<t(end)&&Location(j)>X(1)&&Location(j)<X(end)                   
                    Measure_countArray(floor((Time(j)-t_i)/dt)+1,floor((Location(j)-X_i)/dX)+1)=Measure_countArray(floor((Time(j)-t_i)/dt)+1,floor((Location(j)-X_i)/dX)+1)+1;
                    TotalDensEst(floor((Time(j)-t_i)/dt)+1,floor((Location(j)-X_i)/dX)+1)= TotalDensEst(floor((Time(j)-t_i)/dt)+1,floor((Location(j)-X_i)/dX)+1)+((1/A)*((Velocity(j)+((-1/3)*Acceleration(j)))));
                    TotalPertEst(floor((Time(j)-t_i)/dt)+1,floor((Location(j)-X_i)/dX)+1)=TotalPertEst(floor((Time(j)-t_i)/dt)+1,floor((Location(j)-X_i)/dX)+1)-((-1/3)*(Acceleration(j))/(B));
                    TotalVel(floor((Time(j)-t_i)/dt)+1,floor((Location(j)-X_i)/dX)+1)=TotalVel(floor((Time(j)-t_i)/dt)+1,floor((Location(j)-X_i)/dX)+1)+Velocity(j);
                end 
            end 
            % loop for number of vehicle, 
            for J=1:size(Time,1)
                if Time(J)>t(1)&&Time(J)<t(end)&&Location(J)>X(1)&&Location(J)<X(end)&&J==size(Time,1)
                    VehNum1(floor((Time(size(Time,1))-t_i)/dt)+1,floor((Location(size(Time,1))-X_i)/dX)+1)=VehNum1(floor((Time(size(Time,1))-t_i)/dt)+1,floor((Location(size(Time,1))-X_i)/dX)+1)+1;
                end
                
                if Time(J)>t(1)&&Time(J)<t(end)&&Location(J)>X(1)&&Location(J)<X(end)&&((floor((Time(J)-t_i)/dt)+1)~=(floor((Time(J+1)-t_i)/dt)+1)&&(floor((Location(J)-X_i)/dX)+1)==(floor((Location(J+1)-X_i)/dX)+1))
                    VehNum1(floor((Time(J)-t_i)/dt)+1,floor((Location(J)-X_i)/dX)+1)=VehNum1(floor((Time(J)-t_i)/dt)+1,floor((Location(J)-X_i)/dX)+1)+1;
                end
         
                if Time(J)>t(1)&&Time(J)<t(end)&&Location(J)>X(1)&&Location(J)<X(end)&&((floor((Location(J)-X_i)/dX)+1)~=(floor((Location(J+1)-X_i)/dX)+1))
                    VehNum2(floor((Time(J)-t_i)/dt)+1,floor((Location(J)-X_i)/dX)+1)=(VehNum2(floor((Time(J)-t_i)/dt)+1,floor((Location(J)-X_i)/dX)+1)+1);
                end    
            end
            
            VehNum=VehNum1+(((max(max(VehNum2)))+(min(min(VehNum2)))-VehNum2));
        end
        index_begin=i+1; % change index for starting point of next vehicle
    end

end

%calculate the average velocity and estimated perturbation in each cell
AverVel=TotalVel./Measure_countArray;
PertEst=(TotalPertEst)./Measure_countArray;

%calculation associated and estimated density and velocity of each cell
DensEst=jamDens-(TotalDensEst./Measure_countArray);
VelEst=A*(TotalDensEst./Measure_countArray)+B*(PertEst./DensEst).*(TotalDensEst./Measure_countArray);


%calculating the vehicle density in each cell, the density=number of
%vehicle/the distance
VehDens=VehNum/X_length;

%calculate the vehicle flow in each cell
VehFlow=VehDens.*AverVel;
VehFlowEst=DensEst.*VelEst;

%plot the fundamental diagram of density-flow relationship
figure
plot(VehDens,VehFlow,'r*');
title('Flow-Density 4:00pm-4:15pm');
xlabel('Vehicle Density (veh/foot)');
ylabel('Vehicle Flow (veh/s)');

%plot the density profile pattern
figure
subplot(3,1,1)
image(VehDens','CDataMapping','scale');
colorbar;
title('ground-truth density');
xlabel('Time (second)');
ylabel('Location (foot)');
caxis([0.04,0.27]);

subplot(3,1,2)
image(DensEst','CDataMapping','scale');
colorbar;
title('PTM estimation density');
xlabel('Time (second)');
ylabel('Location (foot)');
caxis([0.04,0.27]);


% starting the Agent-based Modeling process
% loop the TIME t for agents

N=18;%jamDens*X_length;

%creating an zeros matrix for estimated vehicle number
ABMVehNum=zeros(size(VehNum));
ABMVehNum(1,:)=VehNum(1,:);
ABMVehNum(:,end)=VehNum(:,end);
ABMVehNum(:,1)=VehNum(:,1);


[n_row,n_col]=size(VehNum);

for t=1:n_row-1 % iteration as row by row due to the CTMinFlow(t,i+1)
     for i=1:n_col-2
         % If former cell and later cell all have space, then the number of 
         % vehicel in this cell is equal to the original number plus inflow number and minus outflow number, for details, see report. 
         if ABMVehNum(t,i+1)<N&&ABMVehNum(t,i+2)<N
            ABMVehNum(t+1,i+1)=ABMVehNum(t,i+1)+(N-ABMVehNum(t,i+1))-(N-ABMVehNum(t,i+2));
         end
         % the following 3 scenarios all indicate that: if any of the former or later cell is full, the number of cell remains still.   
         if ABMVehNum(t,i+1)<N&&ABMVehNum(t,i+2)>N
            ABMVehNum(t+1,i+1)=ABMVehNum(t,i+1);
         end
         if ABMVehNum(t,i+1)>N&&ABMVehNum(t,i+2)<N
            ABMVehNum(t+1,i+1)=ABMVehNum(t,i+1);
         end
         if ABMVehNum(t,i+1)>N&&ABMVehNum(t,i+2)>N
            ABMVehNum(t+1,i+1)=ABMVehNum(t,i+1);
            
         end
         
     end   
end

% calculating the estimated density of ABM
ABMVehDens=ABMVehNum/X_length;

%plot the density profile pattern and comparing with former two plots
subplot(3,1,3)
image(ABMVehDens','CDataMapping','scale');
colorbar;
title('Agent-based estimation density');
xlabel('Time (second)');
ylabel('Location (foot)');
caxis([0.04,0.27]);

%calculating the relative error of both models
Dif=abs(DensEst-VehDens);
ratio=Dif./VehDens;
sum_1=sum(sum(ratio));
percentrage_PTM=sum_1/(size(VehDens,1)*size(VehDens,2));

Dif2=abs(ABMVehDens-VehDens);
ratio2=Dif2./VehDens;
sum_2=sum(sum(ratio2));
percentrage_ABM=sum_2/((size(VehDens,1)-2)*(size(VehDens,2)-2));


% density analysis according time series
ground=VehDens';
ground_average=mean(ground);

PTM=DensEst';
PTM_average=mean(PTM);

ABM=ABMVehDens';
ABM_average=mean(ABM);

figure
plot(ground_average,'b-');
hold on
plot(PTM_average,'r-');
plot(ABM_average,'k-');
title('Mean density on study area');
xlabel('Time (second)');
ylabel('Density (veh/foot)');
grid on
legend('Ground-true','Agent-based estimation');





