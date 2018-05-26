%function pub_server(varargin)
function pub_server
addpath(genpath('C:\Users\TeleOperation\Downloads\matlab-zmq-master'))

%Solved sign problem still zero and scaling termed removed
%learning task parameters along with accelerations for regression by using mixture of mixtures.
%Also adapting covariance

% clear all;
% close all;
% 
% 
% addpath('C:\Users\TeleOperation\Downloads\TeleoperationData\Experiment1/Subject1');
% addpath('C:\Users\TeleOperation\Downloads\TeleoperationData\Experiment1/Subject2');
% addpath('C:\Users\TeleOperation\Downloads\TeleoperationData\Experiment1/Subject3');
% addpath('C:\Users\TeleOperation\Downloads\TeleoperationData\Experiment1/Subject4');
% addpath('C:\Users\TeleOperation\Downloads\TeleoperationData\Experiment1/Subject5');
% 
% 
% 
% %addpath('/home/affan/Desktop/Papers to read/Sylvain/DMP-learned-by-GMR-v1.0/TeleoperationData/Experiment1/Subject1/');
% %addpath('/home/affan/Desktop/Papers to read/Sylvain/DMP-learned-by-GMR-v1.0/TeleoperationData/Experiment1/Subject2/');
% %addpath('/home/affan/Desktop/Papers to read/Sylvain/DMP-learned-by-GMR-v1.0/TeleoperationData/Experiment1/Subject3/');
% %addpath('/home/affan/Desktop/Papers to read/Sylvain/DMP-learned-by-GMR-v1.0/TeleoperationData/Experiment1/Subject4/');
% %addpath('/home/affan/Desktop/Papers to read/Sylvain/DMP-learned-by-GMR-v1.0/TeleoperationData/Experiment1/Subject5/');
% % addpath('/home/affan/Desktop/Papers to read/Sylvain/DMP-learned-by-GMR-v1.0/ExperimentalDataTeleoperation/Case 02 (Multiple master single slave with predefined authority)//');
% % addpath('/home/affan/Desktop/Papers to read/Sylvain/DMP-learned-by-GMR-v1.0/ExperimentalDataTeleoperation/Case 03 (Multiple master single slave with automatic authority decomposition)//');
% % addpath('/home/affan/Desktop/Papers to read/Sylvain/DMP-learned-by-GMR-v1.0/ExperimentalDataTeleoperation/new_demonstrations/');
% 
% 
% 
% LoadModel = 0; 
% time_duration = 5;
% tao = 1;%1/time_duration;
% Omega = 1/13000;
% wn = 200;
% zeta=0.707;
% kv = 2*zeta*wn; %Damping gain 
% kp = (wn^2)/kv; %Initial stiffness gain
% 
% % currF0 = ((xT-rGMR(nb ).currPos)*kP - rGMR(nb).currVel)*kV
% 
% components = 20; %20 gaussians   
% 
% load('experiment1_subject3_trial1.mat'); 
% load('experiment1_subject3_trial2.mat');
% experiment1_subject3_trial2 = experiment1_subject4_trial1;
% load('experiment1_subject3_trial3.mat');
% load('experiment1_subject3_trial4.mat');
% 
% 
% figure;
% hold on;
% experiment1_subject3_trial1 = experiment1_subject3_trial1(:,2500:19000); 
% plot3(experiment1_subject3_trial1(1,:),experiment1_subject3_trial1(2,:),experiment1_subject3_trial1(3,:));
% figure;
% hold on;
% experiment1_subject3_trial2 = experiment1_subject3_trial2(:,700:end);
% plot3(experiment1_subject3_trial2(1,:),experiment1_subject3_trial2(2,:),experiment1_subject3_trial2(3,:));
% figure;
% hold on;
% experiment1_subject3_trial3 = experiment1_subject3_trial3(:,1:13800);
% plot3(experiment1_subject3_trial3(1,:),experiment1_subject3_trial3(2,:),experiment1_subject3_trial3(3,:));
% figure;
% hold on;
% experiment1_subject3_trial4 = experiment1_subject3_trial4(:,1:13500);
% plot3(experiment1_subject3_trial4(1,:),experiment1_subject3_trial4(2,:),experiment1_subject3_trial4(3,:));
% 
% 
% 
% 
% 
% 
% yy = spline(linspace(0,1,length(experiment1_subject3_trial1)),experiment1_subject3_trial1,linspace(0,1,1300)); %sampled to get 1300 data points for an input time duration of 0 to 1
% yy1=yy'; % for formatting
% yy1_v = (yy1(1:end-1,:)-yy1(2:end,:))/(time_duration/1300); % numerical differentation to get velocity and acceleration, time duration in between two datapoints 
% yy1_a = (yy1_v(1:end-1,:)-yy1_v(2:end,:))/(time_duration/1300);
% 
% yy = spline(linspace(0,1,length(experiment1_subject3_trial2)),experiment1_subject3_trial2,linspace(0,1,1300));
% yy2=yy';
% yy2_v = (yy2(1:end-1,:)-yy2(2:end,:))/(time_duration/1300);
% yy2_a = (yy2_v(1:end-1,:)-yy2_v(2:end,:))/(time_duration/1300);
% 
% 
% yy = spline(linspace(0,1,length(experiment1_subject3_trial3)),experiment1_subject3_trial3,linspace(0,1,1300));
% yy4=yy';
% % yy4=yy4(1:7000,:);
% yy4_v = (yy4(1:end-1,:)-yy4(2:end,:))/(time_duration/1300);
% yy4_a = (yy4_v(1:end-1,:)-yy4_v(2:end,:))/(time_duration/1300);
% 
% 
% yy = spline(linspace(0,1,length(experiment1_subject3_trial4)),experiment1_subject3_trial4,linspace(0,1,1300));
% yy5=yy';
% % yy5=yy5(4000:end,:);
% yy5_v = (yy5(1:end-1,:)-yy5(2:end,:))/(time_duration/1300);
% yy5_a = (yy5_v(1:end-1,:)-yy5_v(2:end,:))/(time_duration/1300);
% 
% 
% g = mean([yy1; yy2; yy4; yy5],1); %for rythmic motion, goal is the mean of the data points in the trajectories 
% yy1r = (yy1_a - tao*kv*(kp*( repmat(g,1298,1) - yy1(1:end-2,:) ) - yy1_v(1:end-1,:) ))./(tao*(repmat(g,1298,1)-repmat(yy1(1,:),length(yy1_a),1))); % F(s) values to be encoded into GMM 
% yy2r = (yy2_a - tao*kv*(kp*( repmat(g,1298,1) - yy2(1:end-2,:) ) - yy2_v(1:end-1,:) ))./(tao*(repmat(g,1298,1)-repmat(yy1(1,:),length(yy2_a),1)));
% % yy3r = (yy3_a - tao*kv*(kp*( repmat(g,length(yy3_a),1) - yy3(1:end-2,:) ) - yy3_v(1:end-1,:) ))./(tao*(repmat(g,length(yy3_a),1)-repmat(yy3(1,:),length(yy3_a),1)));
% yy4r = (yy4_a - tao*kv*(kp*( repmat(g,length(yy4_a),1) - yy4(1:end-2,:) ) - yy4_v(1:end-1,:) ))./(tao*(repmat(g,length(yy4_a),1)-repmat(yy1(1,:),length(yy4_a),1)));
% yy5r = (yy5_a - tao*kv*(kp*( repmat(g,length(yy5_a),1) - yy5(1:end-2,:) ) - yy5_v(1:end-1,:) ))./(tao*(repmat(g,length(yy5_a),1)-repmat(yy1(1,:),length(yy5_a),1)));
% 
% 
% % plot3(Master1_Position_Trial3(1:13284,1),Master1_Position_Trial3(1:13284,2),Master1_Position_Trial3(1:13284,3));
% % hold on;
% % plot3(yy3(:,1),yy3(:,2),yy3(:,3));
% 
% % plot3(Master1_Position_Trial2(1:13000,1),Master1_Position_Trial2(1:13000,2),Master1_Position_Trial2(1:13000,3));
% % plot3(Master1_Position_Trial3(1:13000,1),Master1_Position_Trial3(1:13000,2),Master1_Position_Trial3(1:13000,3));
% % plot3(Master1_Position_Trial4(1:12560,1),Master1_Position_Trial4(1:12560,2),Master1_Position_Trial4(1:12560,3));
% % plot3(Master1_Position_Trial5(1:13400,1),Master1_Position_Trial5(1:13400,2),Master1_Position_Trial5(1:13400,3));
%  
% % figure;
% % hold on;
% % for i=1:10000
% % %     plot3(Master1_Position_Trial1(i,1),Master1_Position_Trial1(i,2),Master1_Position_Trial1(i,3),'or');
% %     plot3(Master1_Position_Trial2(i,1),Master1_Position_Trial2(i,2),Master1_Position_Trial2(i,3),'oy');
% %     plot3(Master1_Position_Trial3(i,1),Master1_Position_Trial3(i,2),Master1_Position_Trial3(i,3),'om');
% %     plot3(Master1_Position_Trial4(i,1),Master1_Position_Trial4(i,2),Master1_Position_Trial4(i,3),'oc');
% %     plot3(Master1_Position_Trial5(i,1),Master1_Position_Trial5(i,2),Master1_Position_Trial5(i,3),'og');   
% %     drawnow;
% % end
% 
% % yy5r=yy5r(5000:end,:);
% % yy4r=yy4r(1:6000,:);
% 
% % PositionData = [Master1_Position_Trial1;Master1_Position_Trial2;Master1_Position_Trial3;Master1_Position_Trial4;Master1_Position_Trial5];
% % PositionData = [[cumsum(0.1*ones(length(Master1_Position_Trial2),1)) Master1_Position_Trial2];[cumsum(0.1*ones(length(Master1_Position_Trial3),1)) Master1_Position_Trial3];[cumsum(0.1*ones(length(Master1_Position_Trial4),1)) Master1_Position_Trial4];[cumsum(0.1*ones(length(Master1_Position_Trial5),1)) Master1_Position_Trial5]];
% 
% % PositionData = [[linspace(0,(2*pi/13000)*12998,12998)' yy1r];[linspace(0,(2*pi/13000)*12998,12998)' yy2r];[linspace(0,(2*pi/13000)*12998,12998)' yy3r];[linspace(0,(2*pi/13000)*12998,12998)' yy4r];[linspace(0,(2*pi/13000)*12998,12998)' yy5r]];
% PositionData = [[linspace(0,(2*pi/1300)*1298,1298)' yy1(1:end-2,:) yy1r]]; % after numerical differentation 2 samples less, first trajectory reference, , yy1 is trajectory, yy1r is F(s)
% % PositionData = PositionData(6000:end,:);
% % PositionData = PositionData(1:4:end,:);
% % Position_Diff_Data_Miss = [[ yy4r];[ yy5r]]';
% Position_Diff_Data_Miss = [ [linspace(0,((2*pi/1300))*12998,length(yy2r))' yy2(1:end-2,:) yy2r];[linspace(0,((2*pi/1300))*12998,length(yy4r))' yy4(1:end-2,:) yy4r]; [linspace(0,((2*pi/1300))*12998,length(yy5r))' yy5(1:end-2,:) yy5r]]';
% % Position_Diff_Data_Miss = [PositionData' Position_Diff_Data_Miss ];
% % Position_Diff_Data_Miss = [[linspace(0,((2*pi/13000))*12998,12998)' yy5r]]';
% % Position_Diff_Data_Miss = Position_Diff_Data_Miss(:,1:4:end);
% % Position_Diff_Data_Miss = [  [linspace(0,(2*pi/13000)*12998,12998)' yy3r]]';
% % PositionData = [[linspace(0,(2*pi/13000)*12998,12998)' yy2r]];
% 
% % Diff_Data = [PositionData(2:end,:)-PositionData(1:end-1,:)];
% 
% % Position_Diff_Data = [PositionData(1:end-1,:) Diff_Data]';
% Position_Diff_Data = [PositionData]';
% 
% % Priors = ones(1,components)/components;
% % ind = ceil(rand(components,1) * size(Position_Diff_Data,2));
% % Mu = Position_Diff_Data(:,ind);
% 
% [Priors, Mu, Sigma] = EM_init_regularTiming(Position_Diff_Data, components); % function: before encoding into GMM, we first initialize the GMM paramters. How? 
% % Sigma = Sigma*10;
% % Sigma = repmat(cov(Position_Diff_Data')/5,1,1,components);
% 
% % [Priors, Mu, Sigma] = EM_boundingCovTele(Position_Diff_Data, Priors, Mu, Sigma);
% % [Priors, Mu, Sigma] = EM_Pol(Position_Diff_Data, 1 , Priors, Mu, Sigma);
% figure
% % plot(yy2r(:,1))
% plot(PositionData(:,1),PositionData(:,5),'*')
% hold on
% plot(Position_Diff_Data_Miss(1,:),Position_Diff_Data_Miss(5,:),'*')
% Miss_Dim = [1]; %missing dimension, phase variable to be estimated s, indexed 
% Obs_Dim = [2 3]; % not to be changed, F(s) and x, indexed
% 
% % Position_Diff_Data = Position_Diff_Data(:,1:50);
% tic
% if(LoadModel)
%     load('DMP_GMM_ModelS2.mat');
% else
%     parfor i=1:3
%         [Priors2{i}, Mu2{i}, Sigma2{i},~, Position_Diff_Data_Miss2{i}] = EM_Pol_MissLocal(Position_Diff_Data([1 (1+i) (4+i)],:),Position_Diff_Data_Miss([1 (1+i) (4+i)],:), 1 , Priors, Mu([1 (1+i) (4+i)],:), Sigma([1 (1+i) (4+i)],[1 (1+i) (4+i)],:), Miss_Dim, Obs_Dim); 
%     end
% %     for i=1:3
% % %        Position_Diff_Data_Miss([(4+i)],:) = Position_Diff_Data_Miss2{i}(end,:); 
% %        [Position_Diff_Data_Miss([(4+i)],:), ~] = GMR_Polar(Priors2{i}, Mu2{i}, Sigma2{i}, Position_Diff_Data_Miss(1,:), 1, [3],1);
% %     end
% end
% toc
% % figure;
% % hold on;
% % for abc=1:size(Mu,2)
% %     plot_gaussian_ellipsoid(Mu2{1}([2 3 4],abc),Sigma2{1}([2 3 4],[2 3 4],abc),1);
% % end 
% 
% % figure;
% % hold on;
% % abc=0:0.1:2*pi;
% % for abc2=2:length(abc)
% %     [val ind] = find(abc(abc2-1)<=Mu(1,:) & Mu(1,:)<=abc(abc2))
% %     for i=ind
% %         plot_gaussian_ellipsoid(Mu([2 3 4],i),Sigma([2 3 4],[2 3 4],i),1);
% %     end
% % end 
% 
% % figure
% % plot(Master1_Position_Trial1(:,1));
% % hold on;
% % plot(Master1_Position_Trial2(:,2));
% % plot(Master1_Position_Trial3(:,3));
% % plot(Master1_Position_Trial4(:,4));
% % plot(Master1_Position_Trial5(:,5));
% 
% % in=[1 2 3];
% % out=[4 5 6];
% in=[1]; %input dimension of GMM
% out=[2 3 4]; %output dimension of GMM
% 
% % Start_Pos = Master1_Position_Trial2(1,:)';
% Clock = 0.0;
% % Traj=[Start_Pos];
% 
% %%model learned above, now we execute the motion
% 
% %predict clock value for random starting point from the trajectories
% RandNum = 3;
% % RandNum = 2;
% % RandomStart = 1;% randi(length(yy3));
% if(RandNum==1)
%     Current_Point = yy1(1,:)';
% %     Current_Point = yy1(randi(length(yy1)),:)';
%     Current_Velocity = (yy1(2,:)-yy1(1,:))'/time_duration;    
% elseif(RandNum==2)
%     Current_Point = yy4(1,:)';
% %     Current_Point = yy4(randi(length(yy4)),:)';
%     Current_Velocity = (yy4(2,:)-yy4(1,:))'/time_duration;
% else
%     Current_Point = yy5(1,:)';
% %     Current_Point = yy5(randi(length(yy5)),:)';
%     Current_Velocity = (yy5(2,:)-yy5(1,:))'/time_duration;
% end
% Start_Point = Current_Point;
% Current_Velocity = Current_Velocity*0;
% % plot3(Start_Point(1),Start_Point(2),Start_Point(3),'*r','MarkerSize',20);
% 
% % PredictedForcingTerm = (Start_Point' - tao*kv*(kp*( g - Start_Point' ) - Current_Velocity' ))./(tao*(g-Start_Point'));
% % [Clock, Sigma_y] = GMR_Polar_max(Priors2{1}, Mu2{1}, Sigma2{1}, Current_Point(1), [2], 1,1);
% 
% Traj=[];
% Traj_Length = 3000;
% 
% Samples = linspace(0,2*pi,300); %clock samples generated by us
% % Probabilities = zeros(3,size(Samples,2));
% 
% % for i=1:length(Current_Point)
% %     Evaluation_Samples = [Samples;repmat(Current_Point(i),1,size(Samples,2))];
% %     for j=1:length(Priors2{i})
% %         Temp = Priors2{i}(j)*gaussPDF([Samples; repmat(Current_Point(i),1,size(Samples,2))],Mu2{i}([1 2],j),Sigma2{i}([1 2],[1 2],j))';
% %         Probabilities(i,:) = Probabilities(i,:) + log(Temp);
% %     end
% % end
% % Probabilities = sum(Probabilities,1);
% % [val ind] = max(Probabilities);
% % Clock = Samples(ind);
% for j=1:length(Samples)
%     for i=1:3
%         [y(i,1), Sigma_y] = GMR_Polar(Priors2{i}, Mu2{i}, Sigma2{i}, Samples(j), 1, [3],1); % predicting F(s) to find out the initial clock signal value
%     end
%     Acc2(j) =  sum(abs(tao*kv*(kp*( g' - Current_Point ) - Current_Velocity ) + tao*(g'-yy1(1,:)').*y));
% end
% [val ind] = min(Acc2);
% Clock = Samples(ind);
% Traj = Current_Point;
% ClockRecord  = [];
% ClockRecord = [ClockRecord Clock];
% for j=1:Traj_Length
%     
% %     Start_Pos = Master1_Position_Trial2(j,:)';
% %     for i=1:length(Priors) 
% %       H(i) = Priors(i) * gaussPDF(Clock,Mu(in,i),Sigma(in,in,i));
% %     end
% %     H = H / sum(H);
% %     F=0;
% %         for i=1:length(Priors)
% %           currFTmp = Mu(out,i) + Sigma(out,in,i)*inv(Sigma(in,in,i)) * (Clock-Mu(in,i));
% %           F = F + H(i) * currFTmp;
% %         end
% %     Start_Pos = F;
%     for i=1:3
%         [y(i,1), Sigma_y] = GMR_Polar(Priors2{i}, Mu2{i}, Sigma2{i}, Clock, 1, [3],1); %for every DMP, we are predicting F(s)
% %         [y(i,1), Sigma_y] = GMR_Polar(Priors2{i}, Mu2{i}, Sigma2{i}, [Clock Current_Point(i)']', [1 2]', [3],1);
%     end
%     Acc =  tao*kv*(kp*( g' - Current_Point ) - Current_Velocity ) + tao*(g'-yy1(1,:)').*y; %DMP equation
%     Current_Velocity = Current_Velocity + Acc*(time_duration/1300); %numerical integration to get velcoity
%     Current_Point = Current_Point + Current_Velocity*(time_duration/1300); %for a given clock value, a position of the end effector is determined
%     
%     Traj = [Traj Current_Point];
%     Clock = Clock + (2*pi/1300);
%     ClockRecord = [ClockRecord Clock];
% end
% % Clock;
% Traj
% 
% figure;
% plot3(Traj(1,:),Traj(2,:),Traj(3,:),'Linewidth',3,'color',[0    0.4470    0.7410]);
% hold on;
% xlim([-50   100]);
% ylim([-60    40]);
% zlim([-100    50]);
% set(gca,'fontsize',30)
% axis([-50,100,-45,20,-90,20])
% view(-50, 30);
% plot3(Start_Point(1),Start_Point(2),Start_Point(3),'o','color',[0    0.4470    0.7410],'MarkerSize',20, 'LineWidth',3,'MarkerFaceColor',[0    0.4470    0.7410]);
% % plot3(Start_Point(1),Start_Point(2),Start_Point(3),'*','color',[0    0.4470    0.7410],'MarkerSize',20);
% % plot3(Master1_Position_Trial1(:,1),Master1_Position_Trial1(:,2),Master1_Position_Trial1(:,3));
% % plot3(Master1_Position_Trial2(:,1),Master1_Position_Trial2(:,2),Master1_Position_Trial2(:,3));
% % figure
% plot3(experiment1_subject3_trial1(1,:),experiment1_subject3_trial1(2,:),experiment1_subject3_trial1(3,:));
% plot3(experiment1_subject3_trial2(1,:),experiment1_subject3_trial2(2,:),experiment1_subject3_trial2(3,:));
% plot3(experiment1_subject3_trial3(1,:),experiment1_subject3_trial3(2,:),experiment1_subject3_trial3(3,:));
% plot3(experiment1_subject3_trial4(1,:),experiment1_subject3_trial4(2,:),experiment1_subject3_trial4(3,:));

% [valC indC] = find(ClockRecord>=2*pi & ClockRecord<=4*pi);
% 
% figure;
% plot(ClockRecord(indC)-2*pi,Traj(1,indC),'Linewidth',3,'color',[0    0.4470    0.7410],'LineWidth',2);
% xlabel('Phase signal')
% ylabel('Position along x-axis [mm]');
% axis([0 2*pi -20 100])
% set(gca,'fontsize',45)
% 
% figure;
% plot(ClockRecord(indC)-2*pi,Traj(2,indC),'Linewidth',3,'color',[0    0.4470    0.7410],'LineWidth',2);
% xlabel('Phase signal')
% ylabel('Position along y-axis [mm]');
% axis([0 2*pi -50 20])
% set(gca,'fontsize',45)
% 
% figure;
% plot(ClockRecord(indC)-2*pi,Traj(3,indC),'Linewidth',3,'color',[0    0.4470    0.7410],'LineWidth',2);
% xlabel('Phase signal')
% ylabel('Position along z-axis [mm]');
% axis([0 2*pi -100 20])
% set(gca,'fontsize',45)
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% hold on;
% plot3(experiment1_subject3_trial1(1,:),experiment1_subject3_trial1(2,:),experiment1_subject3_trial1(3,:),'m','Linewidth',3);
% plot3(experiment1_subject3_trial1(1,1),experiment1_subject3_trial1(2,1),experiment1_subject3_trial1(3,1),'om','MarkerSize',20, 'LineWidth',3,'MarkerFaceColor','m');
% plot3(experiment1_subject3_trial1(1,end),experiment1_subject3_trial1(2,end),experiment1_subject3_trial1(3,end),'sm','MarkerSize',20, 'LineWidth',3,'MarkerFaceColor','m');
% set(gca,'fontsize',30)
% axis([-50,100,-45,20,-90,20])
% view(-50, 30);
% 
% figure
% hold on;
% plot3(experiment1_subject3_trial2(1,:),experiment1_subject3_trial2(2,:),experiment1_subject3_trial2(3,:),'c','Linewidth',3);
% plot3(experiment1_subject3_trial2(1,1),experiment1_subject3_trial2(2,1),experiment1_subject3_trial2(3,1),'oc','MarkerSize',20, 'LineWidth',3,'MarkerFaceColor','c');
% plot3(experiment1_subject3_trial2(1,end),experiment1_subject3_trial2(2,end),experiment1_subject3_trial2(3,end),'sc','MarkerSize',20, 'LineWidth',3,'MarkerFaceColor','c');
% set(gca,'fontsize',30)
% axis([-50,100,-45,20,-90,20])
% view(-50, 30);
% 
% figure
% hold on;
% plot3(experiment1_subject3_trial3(1,:),experiment1_subject3_trial3(2,:),experiment1_subject3_trial3(3,:),'r','Linewidth',3);
% plot3(experiment1_subject3_trial3(1,1),experiment1_subject3_trial3(2,1),experiment1_subject3_trial3(3,1),'or','MarkerSize',20, 'LineWidth',3,'MarkerFaceColor','r');
% plot3(experiment1_subject3_trial3(1,end),experiment1_subject3_trial3(2,end),experiment1_subject3_trial3(3,end),'sr','MarkerSize',20, 'LineWidth',3,'MarkerFaceColor','r');
% set(gca,'fontsize',30)
% axis([-50,100,-45,20,-90,20])
% view(-50, 30);
% 
% figure
% hold on;
% plot3(experiment1_subject3_trial4(1,:),experiment1_subject3_trial4(2,:),experiment1_subject3_trial4(3,:),'g','Linewidth',3);
% plot3(experiment1_subject3_trial4(1,1),experiment1_subject3_trial4(2,1),experiment1_subject3_trial4(3,1),'og','MarkerSize',20, 'LineWidth',3,'MarkerFaceColor','g');
% plot3(experiment1_subject3_trial4(1,end),experiment1_subject3_trial4(2,end),experiment1_subject3_trial4(3,end),'sg','MarkerSize',20, 'LineWidth',3,'MarkerFaceColor','g');
% set(gca,'fontsize',30)
% axis([-50,100,-45,20,-90,20])
% view(-50, 30);
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %All plots
% 
% 
% figure
% hold on;
% plot3(experiment1_subject3_trial1(1,:),experiment1_subject3_trial1(2,:),experiment1_subject3_trial1(3,:),'m','Linewidth',3);
% plot3(experiment1_subject3_trial1(1,1),experiment1_subject3_trial1(2,1),experiment1_subject3_trial1(3,1),'om','MarkerSize',30, 'LineWidth',3,'MarkerFaceColor','m');
% plot3(experiment1_subject3_trial1(1,end),experiment1_subject3_trial1(2,end),experiment1_subject3_trial1(3,end),'sm','MarkerSize',30, 'LineWidth',3,'MarkerFaceColor','m');
% 
% plot3(experiment1_subject3_trial2(1,:),experiment1_subject3_trial2(2,:),experiment1_subject3_trial2(3,:),'c','Linewidth',3);
% plot3(experiment1_subject3_trial2(1,1),experiment1_subject3_trial2(2,1),experiment1_subject3_trial2(3,1),'oc','MarkerSize',30, 'LineWidth',3,'MarkerFaceColor','c');
% plot3(experiment1_subject3_trial2(1,end),experiment1_subject3_trial2(2,end),experiment1_subject3_trial2(3,end),'sc','MarkerSize',30, 'LineWidth',3,'MarkerFaceColor','c');
% 
% plot3(experiment1_subject3_trial3(1,:),experiment1_subject3_trial3(2,:),experiment1_subject3_trial3(3,:),'r','Linewidth',3);
% plot3(experiment1_subject3_trial3(1,1),experiment1_subject3_trial3(2,1),experiment1_subject3_trial3(3,1),'or','MarkerSize',30, 'LineWidth',3,'MarkerFaceColor','r');
% plot3(experiment1_subject3_trial3(1,end),experiment1_subject3_trial3(2,end),experiment1_subject3_trial3(3,end),'sr','MarkerSize',30, 'LineWidth',3,'MarkerFaceColor','r');
% 
% plot3(experiment1_subject3_trial4(1,:),experiment1_subject3_trial4(2,:),experiment1_subject3_trial4(3,:),'g','Linewidth',3);
% plot3(experiment1_subject3_trial4(1,1),experiment1_subject3_trial4(2,1),experiment1_subject3_trial4(3,1),'og','MarkerSize',30, 'LineWidth',3,'MarkerFaceColor','g');
% plot3(experiment1_subject3_trial4(1,end),experiment1_subject3_trial4(2,end),experiment1_subject3_trial4(3,end),'sg','MarkerSize',30, 'LineWidth',3,'MarkerFaceColor','g');
% set(gca,'fontsize',30)
% axis([-50,100,-45,20,-90,20])
% view(-50, 30);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %Unsyncronized trajectories
% 
% 
% figure
% hold on;
% plot(linspace(0,2* pi,length(experiment1_subject3_trial1(1,:))),experiment1_subject3_trial1(1,:),'m','LineWidth',2);
% plot(linspace(0,2* pi,length(experiment1_subject3_trial2(1,:))),experiment1_subject3_trial2(1,:),'c','LineWidth',2);
% plot(linspace(0,2* pi,length(experiment1_subject3_trial3(1,:))),experiment1_subject3_trial3(1,:),'r','LineWidth',2);
% plot(linspace(0,2* pi,length(experiment1_subject3_trial4(1,:))),experiment1_subject3_trial4(1,:),'g','LineWidth',2);
% xlabel('Phase signal')
% ylabel('Position along x-axis [mm]');
% axis([0 2*pi -20 100])
% box on;
% set(gca,'fontsize',45)
% 
% 
% figure
% hold on;
% plot(linspace(0,2* pi,length(experiment1_subject3_trial1(1,:))),experiment1_subject3_trial1(2,:),'m','LineWidth',2);
% plot(linspace(0,2* pi,length(experiment1_subject3_trial2(1,:))),experiment1_subject3_trial2(2,:),'c','LineWidth',2);
% plot(linspace(0,2* pi,length(experiment1_subject3_trial3(1,:))),experiment1_subject3_trial3(2,:),'r','LineWidth',2);
% plot(linspace(0,2* pi,length(experiment1_subject3_trial4(1,:))),experiment1_subject3_trial4(2,:),'g','LineWidth',2);
% xlabel('Phase signal')
% ylabel('Position along y-axis [mm]');
% axis([0 2*pi -50 20])
% box on;
% set(gca,'fontsize',45)
% 
% figure
% hold on;
% plot(linspace(0,2* pi,length(experiment1_subject3_trial1(1,:))),experiment1_subject3_trial1(3,:),'m','LineWidth',2);
% plot(linspace(0,2* pi,length(experiment1_subject3_trial2(1,:))),experiment1_subject3_trial2(3,:),'c','LineWidth',2);
% plot(linspace(0,2* pi,length(experiment1_subject3_trial3(1,:))),experiment1_subject3_trial3(3,:),'r','LineWidth',2);
% plot(linspace(0,2* pi,length(experiment1_subject3_trial4(1,:))),experiment1_subject3_trial4(3,:),'g','LineWidth',2);
% xlabel('Phase signal')
% ylabel('Position along z-axis [mm]');
% axis([0 2*pi -100 20])
% box on;
% set(gca,'fontsize',45)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %Syncronized trajectories
% 
% 
% % figure
% % hold on;
% % plot(linspace(0,2* pi,length(PositionData(:,1))),PositionData(:,2),'m','LineWidth',2);
% % Temp = Position_Diff_Data_Miss2{1}(1,1:1298);
% % Temp2 = Position_Diff_Data_Miss2{1}(2,1:1298); 
% % % [val ind ] = sort(Temp)
% % plot(Temp,Temp2,'*c','LineWidth',2);
% % % plot(Position_Diff_Data_Miss(1,1:1298),Position_Diff_Data_Miss(7,1:1298),'*c');
% % 
% % Temp = Position_Diff_Data_Miss2{1}(1,1+1298:2*1298);
% % Temp2 = Position_Diff_Data_Miss2{1}(2,1+1298:2*1298); 
% % % [val ind ] = sort(Temp)
% % plot(Temp,Temp2,'*r','LineWidth',2);
% % 
% % Temp = Position_Diff_Data_Miss2{1}(1,1+1298*2:3*1298);
% % Temp2 = Position_Diff_Data_Miss2{1}(2,1+1298*2:3*1298); 
% % % [val ind ] = sort(Temp)
% % plot(Temp,Temp2,'*g','LineWidth',2);
% % 
% % xlabel('Phase signal')
% % ylabel('Position along x-axis [mm]');
% % axis([0 2*pi -20 100])
% % box on;
% % set(gca,'fontsize',45)
% % 
% % 
% % figure
% % hold on;
% % plot(linspace(0,2* pi,length(PositionData(:,1))),PositionData(:,3),'m','LineWidth',2);
% % Temp = Position_Diff_Data_Miss2{2}(1,1:1298);
% % Temp2 = Position_Diff_Data_Miss2{2}(2,1:1298); 
% % % [val ind ] = sort(Temp)
% % plot(Temp,Temp2,'*c','LineWidth',2);
% % % plot(Position_Diff_Data_Miss(1,1:1298),Position_Diff_Data_Miss(7,1:1298),'*c');
% % 
% % Temp = Position_Diff_Data_Miss2{2}(1,1+1298:2*1298);
% % Temp2 = Position_Diff_Data_Miss2{2}(2,1+1298:2*1298); 
% % % [val ind ] = sort(Temp)
% % plot(Temp,Temp2,'*r','LineWidth',2);
% % 
% % Temp = Position_Diff_Data_Miss2{2}(1,1+1298*2:3*1298);
% % Temp2 = Position_Diff_Data_Miss2{2}(2,1+1298*2:3*1298); 
% % % [val ind ] = sort(Temp)
% % plot(Temp,Temp2,'*g','LineWidth',2);
% % 
% % xlabel('Phase signal')
% % ylabel('Position along x-axis [mm]');
% % axis([0 2*pi -20 100])
% % box on;
% % set(gca,'fontsize',45)
% % 
% % figure
% % hold on;
% % plot(linspace(0,2* pi,length(PositionData(:,1))),PositionData(:,4),'m','LineWidth',2);
% % Temp = Position_Diff_Data_Miss2{3}(1,1:1298);
% % Temp2 = Position_Diff_Data_Miss2{3}(2,1:1298); 
% % % [val ind ] = sort(Temp)
% % plot(Temp,Temp2,'*c','LineWidth',2);
% % % plot(Position_Diff_Data_Miss(1,1:1298),Position_Diff_Data_Miss(7,1:1298),'*c');
% % 
% % Temp = Position_Diff_Data_Miss2{3}(1,1+1298:2*1298);
% % Temp2 = Position_Diff_Data_Miss2{3}(2,1+1298:2*1298); 
% % % [val ind ] = sort(Temp)
% % plot(Temp,Temp2,'*r','LineWidth',2);
% % 
% % Temp = Position_Diff_Data_Miss2{3}(1,1+1298*2:3*1298);
% % Temp2 = Position_Diff_Data_Miss2{3}(2,1+1298*2:3*1298); 
% % % [val ind ] = sort(Temp)
% % plot(Temp,Temp2,'*g','LineWidth',2);
% % 
% % xlabel('Phase signal')
% % ylabel('Position along x-axis [mm]');
% % axis([0 2*pi -20 100])
% % box on;
% % set(gca,'fontsize',45)
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if(~LoadModel)
%     save('DMP_GMM_ModelS2.mat','Priors2','Mu2', 'Sigma2');
% end
% 
% figure
% % plot(yy2r(:,1))
% hold on
% plot(PositionData(:,1),PositionData(:,5),'*')
% plot(Position_Diff_Data_Miss(1,:),Position_Diff_Data_Miss(5,:),'*')
% figure
% plot(yy1r(:,1))
% hold on
% plot(yy2r(:,1))
% plot(yy4r(:,1))
% plot(yy5r(:,1))
% figure
% hold on
% plotGMM(Mu([1,2],:), Sigma([1,2],[1,2],:), [0 .8 0], 1);

%%
    % Super reliable weather information server
    %
    % Example borrowed from
    % http://learning-0mq-with-pyzmq.readthedocs.org/en/latest/pyzmq/patterns/pubsub.html
    %
    % This weather server will broadcast messages consisting of two integers
    % separated by a space, the first one is the CEP (ZIP brazilian equivalent)
    % and the second is the temperature in celsius

    %port1 = 5566;
%     if (nargin > 0)
%         port =  varargin{1};
%     end


    port = 5577;
    
    context = zmq.core.ctx_new();
%     sync = zmq.core.socket(context, 'ZMQ_PULL');
%     zmq.core.bind(sync, sprintf('tcp://*:%d', port1));

    publisher = zmq.core.socket(context, 'ZMQ_PUB');
    zmq.core.bind(publisher, sprintf('tcp://*:%d', port));
    
    %message=sscanf(char(zmq.core.recv(sync)), '%d');
    %fprintf('data received: %d\n', message);

    
    fprintf('Broadcasting learned Trajectory information...\n');
%     while (1)
%         %topic = randi([15198, 15202], 1, 1); % Choose a brazilian CEP code (first 5 digits)
%         data = randi([10, 45]);              % Pick a random temperature
%         message2 = sprintf('%d %d', data);
%         fprintf('%s\n', message2);
%         zmq.core.send(socket, uint8(message2));
%         delay(1);
%     end

 for i=1:3001
     
        Px=sprintf('%0.4f ',Traj(1,i));
        Py=sprintf('%0.4f ',Traj(2,i));
        Pz=sprintf('%0.4f ',Traj(3,i));
   
        data=[Px Py Pz];
   
    
    fprintf('Sending Position (Px Py Pz) %0.4f % 0.4f % 0.4f to C program\n', Traj(1,i), Traj(2,i), Traj(3,i));
    
        zmq.core.send(publisher, uint8(data));  
     
        %topic = randi([15198, 15202], 1, 1); % Choose a brazilian CEP code (first 5 digits)
        %delay(1);
  end

    zmq.core.disconnect(publisher, sprintf('tcp://*:%d',port));
    zmq.core.close(publisher);
    %zmq.core.disconnect(sync, sprintf('tcp://*:%d',port1));
    %zmq.core.close(sync);
    zmq.core.ctx_shutdown(context);
    zmq.core.ctx_term(context);             
            
            
end

function delay(seconds)
    % pause the program
    %
    % Aguments:
    %  seconds - delay time in seconds
    tic;
    while toc < seconds
    end
end
