clc;clear;close all;
%% Persistent Attributes
Persistent_Multiplicity=[];
Persistent_Mean=[];
Persistent_Standard_Deviation=[];
Persistent_Maximum=[];
Persistent_Minimum=[];
Persistent_Laplacian_Graph_Energy=[];
Persistent_Generalized_Mean_Graph_Energy=[];
Persistent_Moment_Second_Order=[];
Persistent_Number_Non_Zero_Eigenvalue=[];
Persistent_Quasi_Wiener_Index=[];
Persistent_Spanning_Tree_Number=[];
%% H1
%H1 ESC
H1_ESC_Persistent_Multiplicity=[];
H1_ESC_Persistent_Mean=[];
H1_ESC_Persistent_Standard_Deviation=[];
H1_ESC_Persistent_Maximum=[];
H1_ESC_Persistent_Minimum=[];
H1_ESC_Persistent_Laplacian_Graph_Energy=[];
H1_ESC_Persistent_Generalized_Mean_Graph_Energy=[];
H1_ESC_Persistent_Moment_Second_Order=[];
H1_ESC_Persistent_Number_Non_Zero_Eigenvalue=[];
H1_ESC_Persistent_Quasi_Wiener_Index=[];
H1_ESC_Persistent_Spanning_Tree_Number=[];
load('H1_ESC_Chromosome_VR_L0_EV.mat');
H1_ESC_Chr_Num=length(H1_ESC_Chromosome_VR_L0_EV);
for H1_ESC_Chr=1:H1_ESC_Chr_Num
    H1_ESC_Eigenvalues=H1_ESC_Chromosome_VR_L0_EV{H1_ESC_Chr};
    H1_ESC_Eigenvalues(:,1)=[];
    H1_ESC_Eigenvalues(101,:)=[];%Remove the eigenvalues when fully connected
    
    H1_ESC_Number_Zero_Eigenvalues=sum(H1_ESC_Eigenvalues==0,2);
    H1_ESC_Persistent_Multiplicity=[H1_ESC_Persistent_Multiplicity, H1_ESC_Number_Zero_Eigenvalues];
    
    H1_ESC_Eigenvalues(H1_ESC_Eigenvalues==0) = NaN;
    
    H1_ESC_Eigenvalues_Mean=nanmean(H1_ESC_Eigenvalues,2);
    H1_ESC_Eigenvalues_Mean(isnan(H1_ESC_Eigenvalues_Mean)) = 0;
    H1_ESC_Persistent_Mean=[H1_ESC_Persistent_Mean,H1_ESC_Eigenvalues_Mean];
    
    H1_ESC_Eigenvalues_Standard_Deviation=nanstd(H1_ESC_Eigenvalues,0,2);
    H1_ESC_Eigenvalues_Standard_Deviation(isnan(H1_ESC_Eigenvalues_Standard_Deviation)) = 0;
    H1_ESC_Persistent_Standard_Deviation=[H1_ESC_Persistent_Standard_Deviation,H1_ESC_Eigenvalues_Standard_Deviation];
    
    H1_ESC_Eigenvalues_Maximum=nanmax(H1_ESC_Eigenvalues,[],2);
    H1_ESC_Eigenvalues_Maximum(isnan(H1_ESC_Eigenvalues_Maximum)) = 0;
    H1_ESC_Persistent_Maximum=[H1_ESC_Persistent_Maximum,H1_ESC_Eigenvalues_Maximum];
    
    H1_ESC_Eigenvalues_Minimum=nanmin(H1_ESC_Eigenvalues,[],2);
    H1_ESC_Eigenvalues_Minimum(isnan(H1_ESC_Eigenvalues_Minimum)) = 0;
    H1_ESC_Persistent_Minimum=[H1_ESC_Persistent_Minimum,H1_ESC_Eigenvalues_Minimum];
    
    H1_ESC_Eigenvalues_Sum=nansum(H1_ESC_Eigenvalues,2);
    H1_ESC_Eigenvalues_Sum(isnan(H1_ESC_Eigenvalues_Sum)) = 0;
    H1_ESC_Persistent_Laplacian_Graph_Energy=[H1_ESC_Persistent_Laplacian_Graph_Energy,H1_ESC_Eigenvalues_Sum];
    
    H1_ESC_Eigenvalues_Sum_Absolute_Deviation=nansum(abs(H1_ESC_Eigenvalues-nanmean(H1_ESC_Eigenvalues,2)),2);
    H1_ESC_Eigenvalues_Sum_Absolute_Deviation(isnan(H1_ESC_Eigenvalues_Sum_Absolute_Deviation)) = 0;
    H1_ESC_Persistent_Generalized_Mean_Graph_Energy=[H1_ESC_Persistent_Generalized_Mean_Graph_Energy,H1_ESC_Eigenvalues_Sum_Absolute_Deviation];
    
    H1_ESC_Eigenvalues_Moment_Second_Order=nansum(H1_ESC_Eigenvalues.^2,2);
    H1_ESC_Eigenvalues_Moment_Second_Order(isnan(H1_ESC_Eigenvalues_Moment_Second_Order)) = 0;
    H1_ESC_Persistent_Moment_Second_Order=[H1_ESC_Persistent_Moment_Second_Order,H1_ESC_Eigenvalues_Moment_Second_Order];
    
    H1_ESC_Number_Non_Zero_Eigenvalues=sum(isnan(H1_ESC_Eigenvalues)==0,2);
    H1_ESC_Persistent_Number_Non_Zero_Eigenvalue=[H1_ESC_Persistent_Number_Non_Zero_Eigenvalue,H1_ESC_Number_Non_Zero_Eigenvalues];
    
    H1_ESC_Eigenvalues_Quasi_Wiener_Index=nansum((H1_ESC_Number_Non_Zero_Eigenvalues+1)./H1_ESC_Eigenvalues,2);
    H1_ESC_Eigenvalues_Quasi_Wiener_Index(isnan(H1_ESC_Eigenvalues_Quasi_Wiener_Index)) = 0;
    H1_ESC_Persistent_Quasi_Wiener_Index=[H1_ESC_Persistent_Quasi_Wiener_Index,H1_ESC_Eigenvalues_Quasi_Wiener_Index];
    
    H1_ESC_Eigenvalues_Spanning_Tree_Number=-log(H1_ESC_Number_Non_Zero_Eigenvalues+1)+nansum(log(H1_ESC_Eigenvalues),2);
    H1_ESC_Eigenvalues_Spanning_Tree_Number(isnan(H1_ESC_Eigenvalues_Spanning_Tree_Number)) = 0;
    H1_ESC_Persistent_Spanning_Tree_Number=[H1_ESC_Persistent_Spanning_Tree_Number,H1_ESC_Eigenvalues_Spanning_Tree_Number];
end
H1_ESC_Average_Persistent_Multiplicity=mean(H1_ESC_Persistent_Multiplicity,2);
Persistent_Multiplicity=[Persistent_Multiplicity,H1_ESC_Average_Persistent_Multiplicity];

H1_ESC_Average_Persistent_Mean=mean(H1_ESC_Persistent_Mean,2);
Persistent_Mean=[Persistent_Mean,H1_ESC_Average_Persistent_Mean];

H1_ESC_Average_Persistent_Standard_Deviation=mean(H1_ESC_Persistent_Standard_Deviation,2);
Persistent_Standard_Deviation=[Persistent_Standard_Deviation,H1_ESC_Average_Persistent_Standard_Deviation];

H1_ESC_Average_Persistent_Maximum=mean(H1_ESC_Persistent_Maximum,2);
Persistent_Maximum=[Persistent_Maximum,H1_ESC_Average_Persistent_Maximum];

H1_ESC_Average_Persistent_Minimum=mean(H1_ESC_Persistent_Minimum,2);
Persistent_Minimum=[Persistent_Minimum,H1_ESC_Average_Persistent_Minimum];

H1_ESC_Average_Persistent_Laplacian_Graph_Energy=mean(H1_ESC_Persistent_Laplacian_Graph_Energy,2);
Persistent_Laplacian_Graph_Energy=[Persistent_Laplacian_Graph_Energy,H1_ESC_Average_Persistent_Laplacian_Graph_Energy];

H1_ESC_Average_Persistent_Generalized_Mean_Graph_Energy=mean(H1_ESC_Persistent_Generalized_Mean_Graph_Energy,2);
Persistent_Generalized_Mean_Graph_Energy=[Persistent_Generalized_Mean_Graph_Energy,H1_ESC_Average_Persistent_Generalized_Mean_Graph_Energy];

H1_ESC_Average_Persistent_Moment_Second_Order=mean(H1_ESC_Persistent_Moment_Second_Order,2);
Persistent_Moment_Second_Order=[Persistent_Moment_Second_Order,H1_ESC_Average_Persistent_Moment_Second_Order];

H1_ESC_Average_Persistent_Number_Non_Zero_Eigenvalue=mean(H1_ESC_Persistent_Number_Non_Zero_Eigenvalue,2);
Persistent_Number_Non_Zero_Eigenvalue=[Persistent_Number_Non_Zero_Eigenvalue,H1_ESC_Average_Persistent_Number_Non_Zero_Eigenvalue];

H1_ESC_Average_Persistent_Quasi_Wiener_Index=mean(H1_ESC_Persistent_Quasi_Wiener_Index,2);
Persistent_Quasi_Wiener_Index=[Persistent_Quasi_Wiener_Index,H1_ESC_Average_Persistent_Quasi_Wiener_Index];

H1_ESC_Average_Persistent_Spanning_Tree_Number=mean(H1_ESC_Persistent_Spanning_Tree_Number,2);
Persistent_Spanning_Tree_Number=[Persistent_Spanning_Tree_Number,H1_ESC_Average_Persistent_Spanning_Tree_Number];

%H1 ME
H1_ME_Persistent_Multiplicity=[];
H1_ME_Persistent_Mean=[];
H1_ME_Persistent_Standard_Deviation=[];
H1_ME_Persistent_Maximum=[];
H1_ME_Persistent_Minimum=[];
H1_ME_Persistent_Laplacian_Graph_Energy=[];
H1_ME_Persistent_Generalized_Mean_Graph_Energy=[];
H1_ME_Persistent_Moment_Second_Order=[];
H1_ME_Persistent_Number_Non_Zero_Eigenvalue=[];
H1_ME_Persistent_Quasi_Wiener_Index=[];
H1_ME_Persistent_Spanning_Tree_Number=[];
load('H1_ME_Chromosome_VR_L0_EV.mat');
H1_ME_Chr_Num=length(H1_ME_Chromosome_VR_L0_EV);
for H1_ME_Chr=1:H1_ME_Chr_Num
    H1_ME_Eigenvalues=H1_ME_Chromosome_VR_L0_EV{H1_ME_Chr};
    H1_ME_Eigenvalues(:,1)=[];
    H1_ME_Eigenvalues(101,:)=[];%Remove the eigenvalues when fully connected
    
    H1_ME_Number_Zero_Eigenvalues=sum(H1_ME_Eigenvalues==0,2);
    H1_ME_Persistent_Multiplicity=[H1_ME_Persistent_Multiplicity, H1_ME_Number_Zero_Eigenvalues];
    
    H1_ME_Eigenvalues(H1_ME_Eigenvalues==0) = NaN;
    
    H1_ME_Eigenvalues_Mean=nanmean(H1_ME_Eigenvalues,2);
    H1_ME_Eigenvalues_Mean(isnan(H1_ME_Eigenvalues_Mean)) = 0;
    H1_ME_Persistent_Mean=[H1_ME_Persistent_Mean,H1_ME_Eigenvalues_Mean];
    
    H1_ME_Eigenvalues_Standard_Deviation=nanstd(H1_ME_Eigenvalues,0,2);
    H1_ME_Eigenvalues_Standard_Deviation(isnan(H1_ME_Eigenvalues_Standard_Deviation)) = 0;
    H1_ME_Persistent_Standard_Deviation=[H1_ME_Persistent_Standard_Deviation,H1_ME_Eigenvalues_Standard_Deviation];
    
    H1_ME_Eigenvalues_Maximum=nanmax(H1_ME_Eigenvalues,[],2);
    H1_ME_Eigenvalues_Maximum(isnan(H1_ME_Eigenvalues_Maximum)) = 0;
    H1_ME_Persistent_Maximum=[H1_ME_Persistent_Maximum,H1_ME_Eigenvalues_Maximum];
    
    H1_ME_Eigenvalues_Minimum=nanmin(H1_ME_Eigenvalues,[],2);
    H1_ME_Eigenvalues_Minimum(isnan(H1_ME_Eigenvalues_Minimum)) = 0;
    H1_ME_Persistent_Minimum=[H1_ME_Persistent_Minimum,H1_ME_Eigenvalues_Minimum];
    
    H1_ME_Eigenvalues_Sum=nansum(H1_ME_Eigenvalues,2);
    H1_ME_Eigenvalues_Sum(isnan(H1_ME_Eigenvalues_Sum)) = 0;
    H1_ME_Persistent_Laplacian_Graph_Energy=[H1_ME_Persistent_Laplacian_Graph_Energy,H1_ME_Eigenvalues_Sum];
    
    H1_ME_Eigenvalues_Sum_Absolute_Deviation=nansum(abs(H1_ME_Eigenvalues-nanmean(H1_ME_Eigenvalues,2)),2);
    H1_ME_Eigenvalues_Sum_Absolute_Deviation(isnan(H1_ME_Eigenvalues_Sum_Absolute_Deviation)) = 0;
    H1_ME_Persistent_Generalized_Mean_Graph_Energy=[H1_ME_Persistent_Generalized_Mean_Graph_Energy,H1_ME_Eigenvalues_Sum_Absolute_Deviation];
    
    H1_ME_Eigenvalues_Moment_Second_Order=nansum(H1_ME_Eigenvalues.^2,2);
    H1_ME_Eigenvalues_Moment_Second_Order(isnan(H1_ME_Eigenvalues_Moment_Second_Order)) = 0;
    H1_ME_Persistent_Moment_Second_Order=[H1_ME_Persistent_Moment_Second_Order,H1_ME_Eigenvalues_Moment_Second_Order];
    
    H1_ME_Number_Non_Zero_Eigenvalues=sum(isnan(H1_ME_Eigenvalues)==0,2);
    H1_ME_Persistent_Number_Non_Zero_Eigenvalue=[H1_ME_Persistent_Number_Non_Zero_Eigenvalue,H1_ME_Number_Non_Zero_Eigenvalues];
    
    H1_ME_Eigenvalues_Quasi_Wiener_Index=nansum((H1_ME_Number_Non_Zero_Eigenvalues+1)./H1_ME_Eigenvalues,2);
    H1_ME_Eigenvalues_Quasi_Wiener_Index(isnan(H1_ME_Eigenvalues_Quasi_Wiener_Index)) = 0;
    H1_ME_Persistent_Quasi_Wiener_Index=[H1_ME_Persistent_Quasi_Wiener_Index,H1_ME_Eigenvalues_Quasi_Wiener_Index];
    
    H1_ME_Eigenvalues_Spanning_Tree_Number=-log(H1_ME_Number_Non_Zero_Eigenvalues+1)+nansum(log(H1_ME_Eigenvalues),2);
    H1_ME_Eigenvalues_Spanning_Tree_Number(isnan(H1_ME_Eigenvalues_Spanning_Tree_Number)) = 0;
    H1_ME_Persistent_Spanning_Tree_Number=[H1_ME_Persistent_Spanning_Tree_Number,H1_ME_Eigenvalues_Spanning_Tree_Number];
end
H1_ME_Average_Persistent_Multiplicity=mean(H1_ME_Persistent_Multiplicity,2);
Persistent_Multiplicity=[Persistent_Multiplicity,H1_ME_Average_Persistent_Multiplicity];

H1_ME_Average_Persistent_Mean=mean(H1_ME_Persistent_Mean,2);
Persistent_Mean=[Persistent_Mean,H1_ME_Average_Persistent_Mean];

H1_ME_Average_Persistent_Standard_Deviation=mean(H1_ME_Persistent_Standard_Deviation,2);
Persistent_Standard_Deviation=[Persistent_Standard_Deviation,H1_ME_Average_Persistent_Standard_Deviation];

H1_ME_Average_Persistent_Maximum=mean(H1_ME_Persistent_Maximum,2);
Persistent_Maximum=[Persistent_Maximum,H1_ME_Average_Persistent_Maximum];

H1_ME_Average_Persistent_Minimum=mean(H1_ME_Persistent_Minimum,2);
Persistent_Minimum=[Persistent_Minimum,H1_ME_Average_Persistent_Minimum];

H1_ME_Average_Persistent_Laplacian_Graph_Energy=mean(H1_ME_Persistent_Laplacian_Graph_Energy,2);
Persistent_Laplacian_Graph_Energy=[Persistent_Laplacian_Graph_Energy,H1_ME_Average_Persistent_Laplacian_Graph_Energy];

H1_ME_Average_Persistent_Generalized_Mean_Graph_Energy=mean(H1_ME_Persistent_Generalized_Mean_Graph_Energy,2);
Persistent_Generalized_Mean_Graph_Energy=[Persistent_Generalized_Mean_Graph_Energy,H1_ME_Average_Persistent_Generalized_Mean_Graph_Energy];

H1_ME_Average_Persistent_Moment_Second_Order=mean(H1_ME_Persistent_Moment_Second_Order,2);
Persistent_Moment_Second_Order=[Persistent_Moment_Second_Order,H1_ME_Average_Persistent_Moment_Second_Order];

H1_ME_Average_Persistent_Number_Non_Zero_Eigenvalue=mean(H1_ME_Persistent_Number_Non_Zero_Eigenvalue,2);
Persistent_Number_Non_Zero_Eigenvalue=[Persistent_Number_Non_Zero_Eigenvalue,H1_ME_Average_Persistent_Number_Non_Zero_Eigenvalue];

H1_ME_Average_Persistent_Quasi_Wiener_Index=mean(H1_ME_Persistent_Quasi_Wiener_Index,2);
Persistent_Quasi_Wiener_Index=[Persistent_Quasi_Wiener_Index,H1_ME_Average_Persistent_Quasi_Wiener_Index];

H1_ME_Average_Persistent_Spanning_Tree_Number=mean(H1_ME_Persistent_Spanning_Tree_Number,2);
Persistent_Spanning_Tree_Number=[Persistent_Spanning_Tree_Number,H1_ME_Average_Persistent_Spanning_Tree_Number];

% H1 MS
H1_MS_Persistent_Multiplicity=[];
H1_MS_Persistent_Mean=[];
H1_MS_Persistent_Standard_Deviation=[];
H1_MS_Persistent_Maximum=[];
H1_MS_Persistent_Minimum=[];
H1_MS_Persistent_Laplacian_Graph_Energy=[];
H1_MS_Persistent_Generalized_Mean_Graph_Energy=[];
H1_MS_Persistent_Moment_Second_Order=[];
H1_MS_Persistent_Number_Non_Zero_Eigenvalue=[];
H1_MS_Persistent_Quasi_Wiener_Index=[];
H1_MS_Persistent_Spanning_Tree_Number=[];
load('H1_MS_Chromosome_VR_L0_EV.mat');
H1_MS_Chr_Num=length(H1_MS_Chromosome_VR_L0_EV);
for H1_MS_Chr=1:H1_MS_Chr_Num
    H1_MS_Eigenvalues=H1_MS_Chromosome_VR_L0_EV{H1_MS_Chr};
    H1_MS_Eigenvalues(:,1)=[];
    H1_MS_Eigenvalues(101,:)=[];%Remove the eigenvalues when fully connected
    
    H1_MS_Number_Zero_Eigenvalues=sum(H1_MS_Eigenvalues==0,2);
    H1_MS_Persistent_Multiplicity=[H1_MS_Persistent_Multiplicity, H1_MS_Number_Zero_Eigenvalues];
    
    H1_MS_Eigenvalues(H1_MS_Eigenvalues==0) = NaN;
    
    H1_MS_Eigenvalues_Mean=nanmean(H1_MS_Eigenvalues,2);
    H1_MS_Eigenvalues_Mean(isnan(H1_MS_Eigenvalues_Mean)) = 0;
    H1_MS_Persistent_Mean=[H1_MS_Persistent_Mean,H1_MS_Eigenvalues_Mean];
    
    H1_MS_Eigenvalues_Standard_Deviation=nanstd(H1_MS_Eigenvalues,0,2);
    H1_MS_Eigenvalues_Standard_Deviation(isnan(H1_MS_Eigenvalues_Standard_Deviation)) = 0;
    H1_MS_Persistent_Standard_Deviation=[H1_MS_Persistent_Standard_Deviation,H1_MS_Eigenvalues_Standard_Deviation];
    
    H1_MS_Eigenvalues_Maximum=nanmax(H1_MS_Eigenvalues,[],2);
    H1_MS_Eigenvalues_Maximum(isnan(H1_MS_Eigenvalues_Maximum)) = 0;
    H1_MS_Persistent_Maximum=[H1_MS_Persistent_Maximum,H1_MS_Eigenvalues_Maximum];
    
    H1_MS_Eigenvalues_Minimum=nanmin(H1_MS_Eigenvalues,[],2);
    H1_MS_Eigenvalues_Minimum(isnan(H1_MS_Eigenvalues_Minimum)) = 0;
    H1_MS_Persistent_Minimum=[H1_MS_Persistent_Minimum,H1_MS_Eigenvalues_Minimum];
    
    H1_MS_Eigenvalues_Sum=nansum(H1_MS_Eigenvalues,2);
    H1_MS_Eigenvalues_Sum(isnan(H1_MS_Eigenvalues_Sum)) = 0;
    H1_MS_Persistent_Laplacian_Graph_Energy=[H1_MS_Persistent_Laplacian_Graph_Energy,H1_MS_Eigenvalues_Sum];
    
    H1_MS_Eigenvalues_Sum_Absolute_Deviation=nansum(abs(H1_MS_Eigenvalues-nanmean(H1_MS_Eigenvalues,2)),2);
    H1_MS_Eigenvalues_Sum_Absolute_Deviation(isnan(H1_MS_Eigenvalues_Sum_Absolute_Deviation)) = 0;
    H1_MS_Persistent_Generalized_Mean_Graph_Energy=[H1_MS_Persistent_Generalized_Mean_Graph_Energy,H1_MS_Eigenvalues_Sum_Absolute_Deviation];
    
    H1_MS_Eigenvalues_Moment_Second_Order=nansum(H1_MS_Eigenvalues.^2,2);
    H1_MS_Eigenvalues_Moment_Second_Order(isnan(H1_MS_Eigenvalues_Moment_Second_Order)) = 0;
    H1_MS_Persistent_Moment_Second_Order=[H1_MS_Persistent_Moment_Second_Order,H1_MS_Eigenvalues_Moment_Second_Order];
    
    H1_MS_Number_Non_Zero_Eigenvalues=sum(isnan(H1_MS_Eigenvalues)==0,2);
    H1_MS_Persistent_Number_Non_Zero_Eigenvalue=[H1_MS_Persistent_Number_Non_Zero_Eigenvalue,H1_MS_Number_Non_Zero_Eigenvalues];
    
    H1_MS_Eigenvalues_Quasi_Wiener_Index=nansum((H1_MS_Number_Non_Zero_Eigenvalues+1)./H1_MS_Eigenvalues,2);
    H1_MS_Eigenvalues_Quasi_Wiener_Index(isnan(H1_MS_Eigenvalues_Quasi_Wiener_Index)) = 0;
    H1_MS_Persistent_Quasi_Wiener_Index=[H1_MS_Persistent_Quasi_Wiener_Index,H1_MS_Eigenvalues_Quasi_Wiener_Index];
    
    H1_MS_Eigenvalues_Spanning_Tree_Number=-log(H1_MS_Number_Non_Zero_Eigenvalues+1)+nansum(log(H1_MS_Eigenvalues),2);
    H1_MS_Eigenvalues_Spanning_Tree_Number(isnan(H1_MS_Eigenvalues_Spanning_Tree_Number)) = 0;
    H1_MS_Persistent_Spanning_Tree_Number=[H1_MS_Persistent_Spanning_Tree_Number,H1_MS_Eigenvalues_Spanning_Tree_Number];
end
H1_MS_Average_Persistent_Multiplicity=mean(H1_MS_Persistent_Multiplicity,2);
Persistent_Multiplicity=[Persistent_Multiplicity,H1_MS_Average_Persistent_Multiplicity];

H1_MS_Average_Persistent_Mean=mean(H1_MS_Persistent_Mean,2);
Persistent_Mean=[Persistent_Mean,H1_MS_Average_Persistent_Mean];

H1_MS_Average_Persistent_Standard_Deviation=mean(H1_MS_Persistent_Standard_Deviation,2);
Persistent_Standard_Deviation=[Persistent_Standard_Deviation,H1_MS_Average_Persistent_Standard_Deviation];

H1_MS_Average_Persistent_Maximum=mean(H1_MS_Persistent_Maximum,2);
Persistent_Maximum=[Persistent_Maximum,H1_MS_Average_Persistent_Maximum];

H1_MS_Average_Persistent_Minimum=mean(H1_MS_Persistent_Minimum,2);
Persistent_Minimum=[Persistent_Minimum,H1_MS_Average_Persistent_Minimum];

H1_MS_Average_Persistent_Laplacian_Graph_Energy=mean(H1_MS_Persistent_Laplacian_Graph_Energy,2);
Persistent_Laplacian_Graph_Energy=[Persistent_Laplacian_Graph_Energy,H1_MS_Average_Persistent_Laplacian_Graph_Energy];

H1_MS_Average_Persistent_Generalized_Mean_Graph_Energy=mean(H1_MS_Persistent_Generalized_Mean_Graph_Energy,2);
Persistent_Generalized_Mean_Graph_Energy=[Persistent_Generalized_Mean_Graph_Energy,H1_MS_Average_Persistent_Generalized_Mean_Graph_Energy];

H1_MS_Average_Persistent_Moment_Second_Order=mean(H1_MS_Persistent_Moment_Second_Order,2);
Persistent_Moment_Second_Order=[Persistent_Moment_Second_Order,H1_MS_Average_Persistent_Moment_Second_Order];

H1_MS_Average_Persistent_Number_Non_Zero_Eigenvalue=mean(H1_MS_Persistent_Number_Non_Zero_Eigenvalue,2);
Persistent_Number_Non_Zero_Eigenvalue=[Persistent_Number_Non_Zero_Eigenvalue,H1_MS_Average_Persistent_Number_Non_Zero_Eigenvalue];

H1_MS_Average_Persistent_Quasi_Wiener_Index=mean(H1_MS_Persistent_Quasi_Wiener_Index,2);
Persistent_Quasi_Wiener_Index=[Persistent_Quasi_Wiener_Index,H1_MS_Average_Persistent_Quasi_Wiener_Index];

H1_MS_Average_Persistent_Spanning_Tree_Number=mean(H1_MS_Persistent_Spanning_Tree_Number,2);
Persistent_Spanning_Tree_Number=[Persistent_Spanning_Tree_Number,H1_MS_Average_Persistent_Spanning_Tree_Number];

% H1 NP
H1_NP_Persistent_Multiplicity=[];
H1_NP_Persistent_Mean=[];
H1_NP_Persistent_Standard_Deviation=[];
H1_NP_Persistent_Maximum=[];
H1_NP_Persistent_Minimum=[];
H1_NP_Persistent_Laplacian_Graph_Energy=[];
H1_NP_Persistent_Generalized_Mean_Graph_Energy=[];
H1_NP_Persistent_Moment_Second_Order=[];
H1_NP_Persistent_Number_Non_Zero_Eigenvalue=[];
H1_NP_Persistent_Quasi_Wiener_Index=[];
H1_NP_Persistent_Spanning_Tree_Number=[];
load('H1_NP_Chromosome_VR_L0_EV.mat');
H1_NP_Chr_Num=length(H1_NP_Chromosome_VR_L0_EV);
for H1_NP_Chr=1:H1_NP_Chr_Num
    H1_NP_Eigenvalues=H1_NP_Chromosome_VR_L0_EV{H1_NP_Chr};
    H1_NP_Eigenvalues(:,1)=[];
    H1_NP_Eigenvalues(101,:)=[];%Remove the eigenvalues when fully connected
    
    H1_NP_Number_Zero_Eigenvalues=sum(H1_NP_Eigenvalues==0,2);
    H1_NP_Persistent_Multiplicity=[H1_NP_Persistent_Multiplicity, H1_NP_Number_Zero_Eigenvalues];
    
    H1_NP_Eigenvalues(H1_NP_Eigenvalues==0) = NaN;
    
    H1_NP_Eigenvalues_Mean=nanmean(H1_NP_Eigenvalues,2);
    H1_NP_Eigenvalues_Mean(isnan(H1_NP_Eigenvalues_Mean)) = 0;
    H1_NP_Persistent_Mean=[H1_NP_Persistent_Mean,H1_NP_Eigenvalues_Mean];
    
    H1_NP_Eigenvalues_Standard_Deviation=nanstd(H1_NP_Eigenvalues,0,2);
    H1_NP_Eigenvalues_Standard_Deviation(isnan(H1_NP_Eigenvalues_Standard_Deviation)) = 0;
    H1_NP_Persistent_Standard_Deviation=[H1_NP_Persistent_Standard_Deviation,H1_NP_Eigenvalues_Standard_Deviation];
    
    H1_NP_Eigenvalues_Maximum=nanmax(H1_NP_Eigenvalues,[],2);
    H1_NP_Eigenvalues_Maximum(isnan(H1_NP_Eigenvalues_Maximum)) = 0;
    H1_NP_Persistent_Maximum=[H1_NP_Persistent_Maximum,H1_NP_Eigenvalues_Maximum];
    
    H1_NP_Eigenvalues_Minimum=nanmin(H1_NP_Eigenvalues,[],2);
    H1_NP_Eigenvalues_Minimum(isnan(H1_NP_Eigenvalues_Minimum)) = 0;
    H1_NP_Persistent_Minimum=[H1_NP_Persistent_Minimum,H1_NP_Eigenvalues_Minimum];
    
    H1_NP_Eigenvalues_Sum=nansum(H1_NP_Eigenvalues,2);
    H1_NP_Eigenvalues_Sum(isnan(H1_NP_Eigenvalues_Sum)) = 0;
    H1_NP_Persistent_Laplacian_Graph_Energy=[H1_NP_Persistent_Laplacian_Graph_Energy,H1_NP_Eigenvalues_Sum];
    
    H1_NP_Eigenvalues_Sum_Absolute_Deviation=nansum(abs(H1_NP_Eigenvalues-nanmean(H1_NP_Eigenvalues,2)),2);
    H1_NP_Eigenvalues_Sum_Absolute_Deviation(isnan(H1_NP_Eigenvalues_Sum_Absolute_Deviation)) = 0;
    H1_NP_Persistent_Generalized_Mean_Graph_Energy=[H1_NP_Persistent_Generalized_Mean_Graph_Energy,H1_NP_Eigenvalues_Sum_Absolute_Deviation];
    
    H1_NP_Eigenvalues_Moment_Second_Order=nansum(H1_NP_Eigenvalues.^2,2);
    H1_NP_Eigenvalues_Moment_Second_Order(isnan(H1_NP_Eigenvalues_Moment_Second_Order)) = 0;
    H1_NP_Persistent_Moment_Second_Order=[H1_NP_Persistent_Moment_Second_Order,H1_NP_Eigenvalues_Moment_Second_Order];
    
    H1_NP_Number_Non_Zero_Eigenvalues=sum(isnan(H1_NP_Eigenvalues)==0,2);
    H1_NP_Persistent_Number_Non_Zero_Eigenvalue=[H1_NP_Persistent_Number_Non_Zero_Eigenvalue,H1_NP_Number_Non_Zero_Eigenvalues];
    
    H1_NP_Eigenvalues_Quasi_Wiener_Index=nansum((H1_NP_Number_Non_Zero_Eigenvalues+1)./H1_NP_Eigenvalues,2);
    H1_NP_Eigenvalues_Quasi_Wiener_Index(isnan(H1_NP_Eigenvalues_Quasi_Wiener_Index)) = 0;
    H1_NP_Persistent_Quasi_Wiener_Index=[H1_NP_Persistent_Quasi_Wiener_Index,H1_NP_Eigenvalues_Quasi_Wiener_Index];
    
    H1_NP_Eigenvalues_Spanning_Tree_Number=-log(H1_NP_Number_Non_Zero_Eigenvalues+1)+nansum(log(H1_NP_Eigenvalues),2);
    H1_NP_Eigenvalues_Spanning_Tree_Number(isnan(H1_NP_Eigenvalues_Spanning_Tree_Number)) = 0;
    H1_NP_Persistent_Spanning_Tree_Number=[H1_NP_Persistent_Spanning_Tree_Number,H1_NP_Eigenvalues_Spanning_Tree_Number];
end
H1_NP_Average_Persistent_Multiplicity=mean(H1_NP_Persistent_Multiplicity,2);
Persistent_Multiplicity=[Persistent_Multiplicity,H1_NP_Average_Persistent_Multiplicity];

H1_NP_Average_Persistent_Mean=mean(H1_NP_Persistent_Mean,2);
Persistent_Mean=[Persistent_Mean,H1_NP_Average_Persistent_Mean];

H1_NP_Average_Persistent_Standard_Deviation=mean(H1_NP_Persistent_Standard_Deviation,2);
Persistent_Standard_Deviation=[Persistent_Standard_Deviation,H1_NP_Average_Persistent_Standard_Deviation];

H1_NP_Average_Persistent_Maximum=mean(H1_NP_Persistent_Maximum,2);
Persistent_Maximum=[Persistent_Maximum,H1_NP_Average_Persistent_Maximum];

H1_NP_Average_Persistent_Minimum=mean(H1_NP_Persistent_Minimum,2);
Persistent_Minimum=[Persistent_Minimum,H1_NP_Average_Persistent_Minimum];

H1_NP_Average_Persistent_Laplacian_Graph_Energy=mean(H1_NP_Persistent_Laplacian_Graph_Energy,2);
Persistent_Laplacian_Graph_Energy=[Persistent_Laplacian_Graph_Energy,H1_NP_Average_Persistent_Laplacian_Graph_Energy];

H1_NP_Average_Persistent_Generalized_Mean_Graph_Energy=mean(H1_NP_Persistent_Generalized_Mean_Graph_Energy,2);
Persistent_Generalized_Mean_Graph_Energy=[Persistent_Generalized_Mean_Graph_Energy,H1_NP_Average_Persistent_Generalized_Mean_Graph_Energy];

H1_NP_Average_Persistent_Moment_Second_Order=mean(H1_NP_Persistent_Moment_Second_Order,2);
Persistent_Moment_Second_Order=[Persistent_Moment_Second_Order,H1_NP_Average_Persistent_Moment_Second_Order];

H1_NP_Average_Persistent_Number_Non_Zero_Eigenvalue=mean(H1_NP_Persistent_Number_Non_Zero_Eigenvalue,2);
Persistent_Number_Non_Zero_Eigenvalue=[Persistent_Number_Non_Zero_Eigenvalue,H1_NP_Average_Persistent_Number_Non_Zero_Eigenvalue];

H1_NP_Average_Persistent_Quasi_Wiener_Index=mean(H1_NP_Persistent_Quasi_Wiener_Index,2);
Persistent_Quasi_Wiener_Index=[Persistent_Quasi_Wiener_Index,H1_NP_Average_Persistent_Quasi_Wiener_Index];

H1_NP_Average_Persistent_Spanning_Tree_Number=mean(H1_NP_Persistent_Spanning_Tree_Number,2);
Persistent_Spanning_Tree_Number=[Persistent_Spanning_Tree_Number,H1_NP_Average_Persistent_Spanning_Tree_Number];

% H1 TB
H1_TB_Persistent_Multiplicity=[];
H1_TB_Persistent_Mean=[];
H1_TB_Persistent_Standard_Deviation=[];
H1_TB_Persistent_Maximum=[];
H1_TB_Persistent_Minimum=[];
H1_TB_Persistent_Laplacian_Graph_Energy=[];
H1_TB_Persistent_Generalized_Mean_Graph_Energy=[];
H1_TB_Persistent_Moment_Second_Order=[];
H1_TB_Persistent_Number_Non_Zero_Eigenvalue=[];
H1_TB_Persistent_Quasi_Wiener_Index=[];
H1_TB_Persistent_Spanning_Tree_Number=[];
load('H1_TB_Chromosome_VR_L0_EV.mat');
H1_TB_Chr_Num=length(H1_TB_Chromosome_VR_L0_EV);
for H1_TB_Chr=1:H1_TB_Chr_Num
    H1_TB_Eigenvalues=H1_TB_Chromosome_VR_L0_EV{H1_TB_Chr};
    H1_TB_Eigenvalues(:,1)=[];
    H1_TB_Eigenvalues(101,:)=[];%Remove the eigenvalues when fully connected
    
    H1_TB_Number_Zero_Eigenvalues=sum(H1_TB_Eigenvalues==0,2);
    H1_TB_Persistent_Multiplicity=[H1_TB_Persistent_Multiplicity, H1_TB_Number_Zero_Eigenvalues];
    
    H1_TB_Eigenvalues(H1_TB_Eigenvalues==0) = NaN;
    
    H1_TB_Eigenvalues_Mean=nanmean(H1_TB_Eigenvalues,2);
    H1_TB_Eigenvalues_Mean(isnan(H1_TB_Eigenvalues_Mean)) = 0;
    H1_TB_Persistent_Mean=[H1_TB_Persistent_Mean,H1_TB_Eigenvalues_Mean];
    
    H1_TB_Eigenvalues_Standard_Deviation=nanstd(H1_TB_Eigenvalues,0,2);
    H1_TB_Eigenvalues_Standard_Deviation(isnan(H1_TB_Eigenvalues_Standard_Deviation)) = 0;
    H1_TB_Persistent_Standard_Deviation=[H1_TB_Persistent_Standard_Deviation,H1_TB_Eigenvalues_Standard_Deviation];
    
    H1_TB_Eigenvalues_Maximum=nanmax(H1_TB_Eigenvalues,[],2);
    H1_TB_Eigenvalues_Maximum(isnan(H1_TB_Eigenvalues_Maximum)) = 0;
    H1_TB_Persistent_Maximum=[H1_TB_Persistent_Maximum,H1_TB_Eigenvalues_Maximum];
    
    H1_TB_Eigenvalues_Minimum=nanmin(H1_TB_Eigenvalues,[],2);
    H1_TB_Eigenvalues_Minimum(isnan(H1_TB_Eigenvalues_Minimum)) = 0;
    H1_TB_Persistent_Minimum=[H1_TB_Persistent_Minimum,H1_TB_Eigenvalues_Minimum];
    
    H1_TB_Eigenvalues_Sum=nansum(H1_TB_Eigenvalues,2);
    H1_TB_Eigenvalues_Sum(isnan(H1_TB_Eigenvalues_Sum)) = 0;
    H1_TB_Persistent_Laplacian_Graph_Energy=[H1_TB_Persistent_Laplacian_Graph_Energy,H1_TB_Eigenvalues_Sum];
    
    H1_TB_Eigenvalues_Sum_Absolute_Deviation=nansum(abs(H1_TB_Eigenvalues-nanmean(H1_TB_Eigenvalues,2)),2);
    H1_TB_Eigenvalues_Sum_Absolute_Deviation(isnan(H1_TB_Eigenvalues_Sum_Absolute_Deviation)) = 0;
    H1_TB_Persistent_Generalized_Mean_Graph_Energy=[H1_TB_Persistent_Generalized_Mean_Graph_Energy,H1_TB_Eigenvalues_Sum_Absolute_Deviation];
    
    H1_TB_Eigenvalues_Moment_Second_Order=nansum(H1_TB_Eigenvalues.^2,2);
    H1_TB_Eigenvalues_Moment_Second_Order(isnan(H1_TB_Eigenvalues_Moment_Second_Order)) = 0;
    H1_TB_Persistent_Moment_Second_Order=[H1_TB_Persistent_Moment_Second_Order,H1_TB_Eigenvalues_Moment_Second_Order];
    
    H1_TB_Number_Non_Zero_Eigenvalues=sum(isnan(H1_TB_Eigenvalues)==0,2);
    H1_TB_Persistent_Number_Non_Zero_Eigenvalue=[H1_TB_Persistent_Number_Non_Zero_Eigenvalue,H1_TB_Number_Non_Zero_Eigenvalues];
    
    H1_TB_Eigenvalues_Quasi_Wiener_Index=nansum((H1_TB_Number_Non_Zero_Eigenvalues+1)./H1_TB_Eigenvalues,2);
    H1_TB_Eigenvalues_Quasi_Wiener_Index(isnan(H1_TB_Eigenvalues_Quasi_Wiener_Index)) = 0;
    H1_TB_Persistent_Quasi_Wiener_Index=[H1_TB_Persistent_Quasi_Wiener_Index,H1_TB_Eigenvalues_Quasi_Wiener_Index];
    
    H1_TB_Eigenvalues_Spanning_Tree_Number=-log(H1_TB_Number_Non_Zero_Eigenvalues+1)+nansum(log(H1_TB_Eigenvalues),2);
    H1_TB_Eigenvalues_Spanning_Tree_Number(isnan(H1_TB_Eigenvalues_Spanning_Tree_Number)) = 0;
    H1_TB_Persistent_Spanning_Tree_Number=[H1_TB_Persistent_Spanning_Tree_Number,H1_TB_Eigenvalues_Spanning_Tree_Number];
end
H1_TB_Average_Persistent_Multiplicity=mean(H1_TB_Persistent_Multiplicity,2);
Persistent_Multiplicity=[Persistent_Multiplicity,H1_TB_Average_Persistent_Multiplicity];

H1_TB_Average_Persistent_Mean=mean(H1_TB_Persistent_Mean,2);
Persistent_Mean=[Persistent_Mean,H1_TB_Average_Persistent_Mean];

H1_TB_Average_Persistent_Standard_Deviation=mean(H1_TB_Persistent_Standard_Deviation,2);
Persistent_Standard_Deviation=[Persistent_Standard_Deviation,H1_TB_Average_Persistent_Standard_Deviation];

H1_TB_Average_Persistent_Maximum=mean(H1_TB_Persistent_Maximum,2);
Persistent_Maximum=[Persistent_Maximum,H1_TB_Average_Persistent_Maximum];

H1_TB_Average_Persistent_Minimum=mean(H1_TB_Persistent_Minimum,2);
Persistent_Minimum=[Persistent_Minimum,H1_TB_Average_Persistent_Minimum];

H1_TB_Average_Persistent_Laplacian_Graph_Energy=mean(H1_TB_Persistent_Laplacian_Graph_Energy,2);
Persistent_Laplacian_Graph_Energy=[Persistent_Laplacian_Graph_Energy,H1_TB_Average_Persistent_Laplacian_Graph_Energy];

H1_TB_Average_Persistent_Generalized_Mean_Graph_Energy=mean(H1_TB_Persistent_Generalized_Mean_Graph_Energy,2);
Persistent_Generalized_Mean_Graph_Energy=[Persistent_Generalized_Mean_Graph_Energy,H1_TB_Average_Persistent_Generalized_Mean_Graph_Energy];

H1_TB_Average_Persistent_Moment_Second_Order=mean(H1_TB_Persistent_Moment_Second_Order,2);
Persistent_Moment_Second_Order=[Persistent_Moment_Second_Order,H1_TB_Average_Persistent_Moment_Second_Order];

H1_TB_Average_Persistent_Number_Non_Zero_Eigenvalue=mean(H1_TB_Persistent_Number_Non_Zero_Eigenvalue,2);
Persistent_Number_Non_Zero_Eigenvalue=[Persistent_Number_Non_Zero_Eigenvalue,H1_TB_Average_Persistent_Number_Non_Zero_Eigenvalue];

H1_TB_Average_Persistent_Quasi_Wiener_Index=mean(H1_TB_Persistent_Quasi_Wiener_Index,2);
Persistent_Quasi_Wiener_Index=[Persistent_Quasi_Wiener_Index,H1_TB_Average_Persistent_Quasi_Wiener_Index];

H1_TB_Average_Persistent_Spanning_Tree_Number=mean(H1_TB_Persistent_Spanning_Tree_Number,2);
Persistent_Spanning_Tree_Number=[Persistent_Spanning_Tree_Number,H1_TB_Average_Persistent_Spanning_Tree_Number];

%% RUES2
% RUES2 CM
RUES2_CM_Persistent_Multiplicity=[];
RUES2_CM_Persistent_Mean=[];
RUES2_CM_Persistent_Standard_Deviation=[];
RUES2_CM_Persistent_Maximum=[];
RUES2_CM_Persistent_Minimum=[];
RUES2_CM_Persistent_Laplacian_Graph_Energy=[];
RUES2_CM_Persistent_Generalized_Mean_Graph_Energy=[];
RUES2_CM_Persistent_Moment_Second_Order=[];
RUES2_CM_Persistent_Number_Non_Zero_Eigenvalue=[];
RUES2_CM_Persistent_Quasi_Wiener_Index=[];
RUES2_CM_Persistent_Spanning_Tree_Number=[];
load('RUES2_CM_Chromosome_VR_L0_EV.mat');
RUES2_CM_Chr_Num=length(RUES2_CM_Chromosome_VR_L0_EV);
for RUES2_CM_Chr=1:RUES2_CM_Chr_Num
    RUES2_CM_Eigenvalues=RUES2_CM_Chromosome_VR_L0_EV{RUES2_CM_Chr};
    RUES2_CM_Eigenvalues(:,1)=[];
    RUES2_CM_Eigenvalues(101,:)=[];%Remove the eigenvalues when fully connected
    
    RUES2_CM_Number_Zero_Eigenvalues=sum(RUES2_CM_Eigenvalues==0,2);
    RUES2_CM_Persistent_Multiplicity=[RUES2_CM_Persistent_Multiplicity, RUES2_CM_Number_Zero_Eigenvalues];
    
    RUES2_CM_Eigenvalues(RUES2_CM_Eigenvalues==0) = NaN;
    
    RUES2_CM_Eigenvalues_Mean=nanmean(RUES2_CM_Eigenvalues,2);
    RUES2_CM_Eigenvalues_Mean(isnan(RUES2_CM_Eigenvalues_Mean)) = 0;
    RUES2_CM_Persistent_Mean=[RUES2_CM_Persistent_Mean,RUES2_CM_Eigenvalues_Mean];
    
    RUES2_CM_Eigenvalues_Standard_Deviation=nanstd(RUES2_CM_Eigenvalues,0,2);
    RUES2_CM_Eigenvalues_Standard_Deviation(isnan(RUES2_CM_Eigenvalues_Standard_Deviation)) = 0;
    RUES2_CM_Persistent_Standard_Deviation=[RUES2_CM_Persistent_Standard_Deviation,RUES2_CM_Eigenvalues_Standard_Deviation];
    
    RUES2_CM_Eigenvalues_Maximum=nanmax(RUES2_CM_Eigenvalues,[],2);
    RUES2_CM_Eigenvalues_Maximum(isnan(RUES2_CM_Eigenvalues_Maximum)) = 0;
    RUES2_CM_Persistent_Maximum=[RUES2_CM_Persistent_Maximum,RUES2_CM_Eigenvalues_Maximum];
    
    RUES2_CM_Eigenvalues_Minimum=nanmin(RUES2_CM_Eigenvalues,[],2);
    RUES2_CM_Eigenvalues_Minimum(isnan(RUES2_CM_Eigenvalues_Minimum)) = 0;
    RUES2_CM_Persistent_Minimum=[RUES2_CM_Persistent_Minimum,RUES2_CM_Eigenvalues_Minimum];
    
    RUES2_CM_Eigenvalues_Sum=nansum(RUES2_CM_Eigenvalues,2);
    RUES2_CM_Eigenvalues_Sum(isnan(RUES2_CM_Eigenvalues_Sum)) = 0;
    RUES2_CM_Persistent_Laplacian_Graph_Energy=[RUES2_CM_Persistent_Laplacian_Graph_Energy,RUES2_CM_Eigenvalues_Sum];
    
    RUES2_CM_Eigenvalues_Sum_Absolute_Deviation=nansum(abs(RUES2_CM_Eigenvalues-nanmean(RUES2_CM_Eigenvalues,2)),2);
    RUES2_CM_Eigenvalues_Sum_Absolute_Deviation(isnan(RUES2_CM_Eigenvalues_Sum_Absolute_Deviation)) = 0;
    RUES2_CM_Persistent_Generalized_Mean_Graph_Energy=[RUES2_CM_Persistent_Generalized_Mean_Graph_Energy,RUES2_CM_Eigenvalues_Sum_Absolute_Deviation];
    
    RUES2_CM_Eigenvalues_Moment_Second_Order=nansum(RUES2_CM_Eigenvalues.^2,2);
    RUES2_CM_Eigenvalues_Moment_Second_Order(isnan(RUES2_CM_Eigenvalues_Moment_Second_Order)) = 0;
    RUES2_CM_Persistent_Moment_Second_Order=[RUES2_CM_Persistent_Moment_Second_Order,RUES2_CM_Eigenvalues_Moment_Second_Order];
    
    RUES2_CM_Number_Non_Zero_Eigenvalues=sum(isnan(RUES2_CM_Eigenvalues)==0,2);
    RUES2_CM_Persistent_Number_Non_Zero_Eigenvalue=[RUES2_CM_Persistent_Number_Non_Zero_Eigenvalue,RUES2_CM_Number_Non_Zero_Eigenvalues];
    
    RUES2_CM_Eigenvalues_Quasi_Wiener_Index=nansum((RUES2_CM_Number_Non_Zero_Eigenvalues+1)./RUES2_CM_Eigenvalues,2);
    RUES2_CM_Eigenvalues_Quasi_Wiener_Index(isnan(RUES2_CM_Eigenvalues_Quasi_Wiener_Index)) = 0;
    RUES2_CM_Persistent_Quasi_Wiener_Index=[RUES2_CM_Persistent_Quasi_Wiener_Index,RUES2_CM_Eigenvalues_Quasi_Wiener_Index];
    
    RUES2_CM_Eigenvalues_Spanning_Tree_Number=-log(RUES2_CM_Number_Non_Zero_Eigenvalues+1)+nansum(log(RUES2_CM_Eigenvalues),2);
    RUES2_CM_Eigenvalues_Spanning_Tree_Number(isnan(RUES2_CM_Eigenvalues_Spanning_Tree_Number)) = 0;
    RUES2_CM_Persistent_Spanning_Tree_Number=[RUES2_CM_Persistent_Spanning_Tree_Number,RUES2_CM_Eigenvalues_Spanning_Tree_Number];
end
RUES2_CM_Average_Persistent_Multiplicity=mean(RUES2_CM_Persistent_Multiplicity,2);
Persistent_Multiplicity=[Persistent_Multiplicity,RUES2_CM_Average_Persistent_Multiplicity];

RUES2_CM_Average_Persistent_Mean=mean(RUES2_CM_Persistent_Mean,2);
Persistent_Mean=[Persistent_Mean,RUES2_CM_Average_Persistent_Mean];

RUES2_CM_Average_Persistent_Standard_Deviation=mean(RUES2_CM_Persistent_Standard_Deviation,2);
Persistent_Standard_Deviation=[Persistent_Standard_Deviation,RUES2_CM_Average_Persistent_Standard_Deviation];

RUES2_CM_Average_Persistent_Maximum=mean(RUES2_CM_Persistent_Maximum,2);
Persistent_Maximum=[Persistent_Maximum,RUES2_CM_Average_Persistent_Maximum];

RUES2_CM_Average_Persistent_Minimum=mean(RUES2_CM_Persistent_Minimum,2);
Persistent_Minimum=[Persistent_Minimum,RUES2_CM_Average_Persistent_Minimum];

RUES2_CM_Average_Persistent_Laplacian_Graph_Energy=mean(RUES2_CM_Persistent_Laplacian_Graph_Energy,2);
Persistent_Laplacian_Graph_Energy=[Persistent_Laplacian_Graph_Energy,RUES2_CM_Average_Persistent_Laplacian_Graph_Energy];

RUES2_CM_Average_Persistent_Generalized_Mean_Graph_Energy=mean(RUES2_CM_Persistent_Generalized_Mean_Graph_Energy,2);
Persistent_Generalized_Mean_Graph_Energy=[Persistent_Generalized_Mean_Graph_Energy,RUES2_CM_Average_Persistent_Generalized_Mean_Graph_Energy];

RUES2_CM_Average_Persistent_Moment_Second_Order=mean(RUES2_CM_Persistent_Moment_Second_Order,2);
Persistent_Moment_Second_Order=[Persistent_Moment_Second_Order,RUES2_CM_Average_Persistent_Moment_Second_Order];

RUES2_CM_Average_Persistent_Number_Non_Zero_Eigenvalue=mean(RUES2_CM_Persistent_Number_Non_Zero_Eigenvalue,2);
Persistent_Number_Non_Zero_Eigenvalue=[Persistent_Number_Non_Zero_Eigenvalue,RUES2_CM_Average_Persistent_Number_Non_Zero_Eigenvalue];

RUES2_CM_Average_Persistent_Quasi_Wiener_Index=mean(RUES2_CM_Persistent_Quasi_Wiener_Index,2);
Persistent_Quasi_Wiener_Index=[Persistent_Quasi_Wiener_Index,RUES2_CM_Average_Persistent_Quasi_Wiener_Index];

RUES2_CM_Average_Persistent_Spanning_Tree_Number=mean(RUES2_CM_Persistent_Spanning_Tree_Number,2);
Persistent_Spanning_Tree_Number=[Persistent_Spanning_Tree_Number,RUES2_CM_Average_Persistent_Spanning_Tree_Number];

% RUES2 CP
RUES2_CP_Persistent_Multiplicity=[];
RUES2_CP_Persistent_Mean=[];
RUES2_CP_Persistent_Standard_Deviation=[];
RUES2_CP_Persistent_Maximum=[];
RUES2_CP_Persistent_Minimum=[];
RUES2_CP_Persistent_Laplacian_Graph_Energy=[];
RUES2_CP_Persistent_Generalized_Mean_Graph_Energy=[];
RUES2_CP_Persistent_Moment_Second_Order=[];
RUES2_CP_Persistent_Number_Non_Zero_Eigenvalue=[];
RUES2_CP_Persistent_Quasi_Wiener_Index=[];
RUES2_CP_Persistent_Spanning_Tree_Number=[];
load('RUES2_CP_Chromosome_VR_L0_EV.mat');
RUES2_CP_Chr_Num=length(RUES2_CP_Chromosome_VR_L0_EV);
for RUES2_CP_Chr=1:RUES2_CP_Chr_Num
    RUES2_CP_Eigenvalues=RUES2_CP_Chromosome_VR_L0_EV{RUES2_CP_Chr};
    RUES2_CP_Eigenvalues(:,1)=[];
    RUES2_CP_Eigenvalues(101,:)=[];%Remove the eigenvalues when fully connected
    
    RUES2_CP_Number_Zero_Eigenvalues=sum(RUES2_CP_Eigenvalues==0,2);
    RUES2_CP_Persistent_Multiplicity=[RUES2_CP_Persistent_Multiplicity, RUES2_CP_Number_Zero_Eigenvalues];
    
    RUES2_CP_Eigenvalues(RUES2_CP_Eigenvalues==0) = NaN;
    
    RUES2_CP_Eigenvalues_Mean=nanmean(RUES2_CP_Eigenvalues,2);
    RUES2_CP_Eigenvalues_Mean(isnan(RUES2_CP_Eigenvalues_Mean)) = 0;
    RUES2_CP_Persistent_Mean=[RUES2_CP_Persistent_Mean,RUES2_CP_Eigenvalues_Mean];
    
    RUES2_CP_Eigenvalues_Standard_Deviation=nanstd(RUES2_CP_Eigenvalues,0,2);
    RUES2_CP_Eigenvalues_Standard_Deviation(isnan(RUES2_CP_Eigenvalues_Standard_Deviation)) = 0;
    RUES2_CP_Persistent_Standard_Deviation=[RUES2_CP_Persistent_Standard_Deviation,RUES2_CP_Eigenvalues_Standard_Deviation];
    
    RUES2_CP_Eigenvalues_Maximum=nanmax(RUES2_CP_Eigenvalues,[],2);
    RUES2_CP_Eigenvalues_Maximum(isnan(RUES2_CP_Eigenvalues_Maximum)) = 0;
    RUES2_CP_Persistent_Maximum=[RUES2_CP_Persistent_Maximum,RUES2_CP_Eigenvalues_Maximum];
    
    RUES2_CP_Eigenvalues_Minimum=nanmin(RUES2_CP_Eigenvalues,[],2);
    RUES2_CP_Eigenvalues_Minimum(isnan(RUES2_CP_Eigenvalues_Minimum)) = 0;
    RUES2_CP_Persistent_Minimum=[RUES2_CP_Persistent_Minimum,RUES2_CP_Eigenvalues_Minimum];
    
    RUES2_CP_Eigenvalues_Sum=nansum(RUES2_CP_Eigenvalues,2);
    RUES2_CP_Eigenvalues_Sum(isnan(RUES2_CP_Eigenvalues_Sum)) = 0;
    RUES2_CP_Persistent_Laplacian_Graph_Energy=[RUES2_CP_Persistent_Laplacian_Graph_Energy,RUES2_CP_Eigenvalues_Sum];
    
    RUES2_CP_Eigenvalues_Sum_Absolute_Deviation=nansum(abs(RUES2_CP_Eigenvalues-nanmean(RUES2_CP_Eigenvalues,2)),2);
    RUES2_CP_Eigenvalues_Sum_Absolute_Deviation(isnan(RUES2_CP_Eigenvalues_Sum_Absolute_Deviation)) = 0;
    RUES2_CP_Persistent_Generalized_Mean_Graph_Energy=[RUES2_CP_Persistent_Generalized_Mean_Graph_Energy,RUES2_CP_Eigenvalues_Sum_Absolute_Deviation];
    
    RUES2_CP_Eigenvalues_Moment_Second_Order=nansum(RUES2_CP_Eigenvalues.^2,2);
    RUES2_CP_Eigenvalues_Moment_Second_Order(isnan(RUES2_CP_Eigenvalues_Moment_Second_Order)) = 0;
    RUES2_CP_Persistent_Moment_Second_Order=[RUES2_CP_Persistent_Moment_Second_Order,RUES2_CP_Eigenvalues_Moment_Second_Order];
    
    RUES2_CP_Number_Non_Zero_Eigenvalues=sum(isnan(RUES2_CP_Eigenvalues)==0,2);
    RUES2_CP_Persistent_Number_Non_Zero_Eigenvalue=[RUES2_CP_Persistent_Number_Non_Zero_Eigenvalue,RUES2_CP_Number_Non_Zero_Eigenvalues];
    
    RUES2_CP_Eigenvalues_Quasi_Wiener_Index=nansum((RUES2_CP_Number_Non_Zero_Eigenvalues+1)./RUES2_CP_Eigenvalues,2);
    RUES2_CP_Eigenvalues_Quasi_Wiener_Index(isnan(RUES2_CP_Eigenvalues_Quasi_Wiener_Index)) = 0;
    RUES2_CP_Persistent_Quasi_Wiener_Index=[RUES2_CP_Persistent_Quasi_Wiener_Index,RUES2_CP_Eigenvalues_Quasi_Wiener_Index];
    
    RUES2_CP_Eigenvalues_Spanning_Tree_Number=-log(RUES2_CP_Number_Non_Zero_Eigenvalues+1)+nansum(log(RUES2_CP_Eigenvalues),2);
    RUES2_CP_Eigenvalues_Spanning_Tree_Number(isnan(RUES2_CP_Eigenvalues_Spanning_Tree_Number)) = 0;
    RUES2_CP_Persistent_Spanning_Tree_Number=[RUES2_CP_Persistent_Spanning_Tree_Number,RUES2_CP_Eigenvalues_Spanning_Tree_Number];
end
RUES2_CP_Average_Persistent_Multiplicity=mean(RUES2_CP_Persistent_Multiplicity,2);
Persistent_Multiplicity=[Persistent_Multiplicity,RUES2_CP_Average_Persistent_Multiplicity];

RUES2_CP_Average_Persistent_Mean=mean(RUES2_CP_Persistent_Mean,2);
Persistent_Mean=[Persistent_Mean,RUES2_CP_Average_Persistent_Mean];

RUES2_CP_Average_Persistent_Standard_Deviation=mean(RUES2_CP_Persistent_Standard_Deviation,2);
Persistent_Standard_Deviation=[Persistent_Standard_Deviation,RUES2_CP_Average_Persistent_Standard_Deviation];

RUES2_CP_Average_Persistent_Maximum=mean(RUES2_CP_Persistent_Maximum,2);
Persistent_Maximum=[Persistent_Maximum,RUES2_CP_Average_Persistent_Maximum];

RUES2_CP_Average_Persistent_Minimum=mean(RUES2_CP_Persistent_Minimum,2);
Persistent_Minimum=[Persistent_Minimum,RUES2_CP_Average_Persistent_Minimum];

RUES2_CP_Average_Persistent_Laplacian_Graph_Energy=mean(RUES2_CP_Persistent_Laplacian_Graph_Energy,2);
Persistent_Laplacian_Graph_Energy=[Persistent_Laplacian_Graph_Energy,RUES2_CP_Average_Persistent_Laplacian_Graph_Energy];

RUES2_CP_Average_Persistent_Generalized_Mean_Graph_Energy=mean(RUES2_CP_Persistent_Generalized_Mean_Graph_Energy,2);
Persistent_Generalized_Mean_Graph_Energy=[Persistent_Generalized_Mean_Graph_Energy,RUES2_CP_Average_Persistent_Generalized_Mean_Graph_Energy];

RUES2_CP_Average_Persistent_Moment_Second_Order=mean(RUES2_CP_Persistent_Moment_Second_Order,2);
Persistent_Moment_Second_Order=[Persistent_Moment_Second_Order,RUES2_CP_Average_Persistent_Moment_Second_Order];

RUES2_CP_Average_Persistent_Number_Non_Zero_Eigenvalue=mean(RUES2_CP_Persistent_Number_Non_Zero_Eigenvalue,2);
Persistent_Number_Non_Zero_Eigenvalue=[Persistent_Number_Non_Zero_Eigenvalue,RUES2_CP_Average_Persistent_Number_Non_Zero_Eigenvalue];

RUES2_CP_Average_Persistent_Quasi_Wiener_Index=mean(RUES2_CP_Persistent_Quasi_Wiener_Index,2);
Persistent_Quasi_Wiener_Index=[Persistent_Quasi_Wiener_Index,RUES2_CP_Average_Persistent_Quasi_Wiener_Index];

RUES2_CP_Average_Persistent_Spanning_Tree_Number=mean(RUES2_CP_Persistent_Spanning_Tree_Number,2);
Persistent_Spanning_Tree_Number=[Persistent_Spanning_Tree_Number,RUES2_CP_Average_Persistent_Spanning_Tree_Number];

% RUES2 ESC
RUES2_ESC_Persistent_Multiplicity=[];
RUES2_ESC_Persistent_Mean=[];
RUES2_ESC_Persistent_Standard_Deviation=[];
RUES2_ESC_Persistent_Maximum=[];
RUES2_ESC_Persistent_Minimum=[];
RUES2_ESC_Persistent_Laplacian_Graph_Energy=[];
RUES2_ESC_Persistent_Generalized_Mean_Graph_Energy=[];
RUES2_ESC_Persistent_Moment_Second_Order=[];
RUES2_ESC_Persistent_Number_Non_Zero_Eigenvalue=[];
RUES2_ESC_Persistent_Quasi_Wiener_Index=[];
RUES2_ESC_Persistent_Spanning_Tree_Number=[];
load('RUES2_ESC_Chromosome_VR_L0_EV.mat');
RUES2_ESC_Chr_Num=length(RUES2_ESC_Chromosome_VR_L0_EV);
for RUES2_ESC_Chr=1:RUES2_ESC_Chr_Num
    RUES2_ESC_Eigenvalues=RUES2_ESC_Chromosome_VR_L0_EV{RUES2_ESC_Chr};
    RUES2_ESC_Eigenvalues(:,1)=[];
    RUES2_ESC_Eigenvalues(101,:)=[];%Remove the eigenvalues when fully connected
    
    RUES2_ESC_Number_Zero_Eigenvalues=sum(RUES2_ESC_Eigenvalues==0,2);
    RUES2_ESC_Persistent_Multiplicity=[RUES2_ESC_Persistent_Multiplicity, RUES2_ESC_Number_Zero_Eigenvalues];
    
    RUES2_ESC_Eigenvalues(RUES2_ESC_Eigenvalues==0) = NaN;
    
    RUES2_ESC_Eigenvalues_Mean=nanmean(RUES2_ESC_Eigenvalues,2);
    RUES2_ESC_Eigenvalues_Mean(isnan(RUES2_ESC_Eigenvalues_Mean)) = 0;
    RUES2_ESC_Persistent_Mean=[RUES2_ESC_Persistent_Mean,RUES2_ESC_Eigenvalues_Mean];
    
    RUES2_ESC_Eigenvalues_Standard_Deviation=nanstd(RUES2_ESC_Eigenvalues,0,2);
    RUES2_ESC_Eigenvalues_Standard_Deviation(isnan(RUES2_ESC_Eigenvalues_Standard_Deviation)) = 0;
    RUES2_ESC_Persistent_Standard_Deviation=[RUES2_ESC_Persistent_Standard_Deviation,RUES2_ESC_Eigenvalues_Standard_Deviation];
    
    RUES2_ESC_Eigenvalues_Maximum=nanmax(RUES2_ESC_Eigenvalues,[],2);
    RUES2_ESC_Eigenvalues_Maximum(isnan(RUES2_ESC_Eigenvalues_Maximum)) = 0;
    RUES2_ESC_Persistent_Maximum=[RUES2_ESC_Persistent_Maximum,RUES2_ESC_Eigenvalues_Maximum];
    
    RUES2_ESC_Eigenvalues_Minimum=nanmin(RUES2_ESC_Eigenvalues,[],2);
    RUES2_ESC_Eigenvalues_Minimum(isnan(RUES2_ESC_Eigenvalues_Minimum)) = 0;
    RUES2_ESC_Persistent_Minimum=[RUES2_ESC_Persistent_Minimum,RUES2_ESC_Eigenvalues_Minimum];
    
    RUES2_ESC_Eigenvalues_Sum=nansum(RUES2_ESC_Eigenvalues,2);
    RUES2_ESC_Eigenvalues_Sum(isnan(RUES2_ESC_Eigenvalues_Sum)) = 0;
    RUES2_ESC_Persistent_Laplacian_Graph_Energy=[RUES2_ESC_Persistent_Laplacian_Graph_Energy,RUES2_ESC_Eigenvalues_Sum];
    
    RUES2_ESC_Eigenvalues_Sum_Absolute_Deviation=nansum(abs(RUES2_ESC_Eigenvalues-nanmean(RUES2_ESC_Eigenvalues,2)),2);
    RUES2_ESC_Eigenvalues_Sum_Absolute_Deviation(isnan(RUES2_ESC_Eigenvalues_Sum_Absolute_Deviation)) = 0;
    RUES2_ESC_Persistent_Generalized_Mean_Graph_Energy=[RUES2_ESC_Persistent_Generalized_Mean_Graph_Energy,RUES2_ESC_Eigenvalues_Sum_Absolute_Deviation];
    
    RUES2_ESC_Eigenvalues_Moment_Second_Order=nansum(RUES2_ESC_Eigenvalues.^2,2);
    RUES2_ESC_Eigenvalues_Moment_Second_Order(isnan(RUES2_ESC_Eigenvalues_Moment_Second_Order)) = 0;
    RUES2_ESC_Persistent_Moment_Second_Order=[RUES2_ESC_Persistent_Moment_Second_Order,RUES2_ESC_Eigenvalues_Moment_Second_Order];
    
    RUES2_ESC_Number_Non_Zero_Eigenvalues=sum(isnan(RUES2_ESC_Eigenvalues)==0,2);
    RUES2_ESC_Persistent_Number_Non_Zero_Eigenvalue=[RUES2_ESC_Persistent_Number_Non_Zero_Eigenvalue,RUES2_ESC_Number_Non_Zero_Eigenvalues];
    
    RUES2_ESC_Eigenvalues_Quasi_Wiener_Index=nansum((RUES2_ESC_Number_Non_Zero_Eigenvalues+1)./RUES2_ESC_Eigenvalues,2);
    RUES2_ESC_Eigenvalues_Quasi_Wiener_Index(isnan(RUES2_ESC_Eigenvalues_Quasi_Wiener_Index)) = 0;
    RUES2_ESC_Persistent_Quasi_Wiener_Index=[RUES2_ESC_Persistent_Quasi_Wiener_Index,RUES2_ESC_Eigenvalues_Quasi_Wiener_Index];
    
    RUES2_ESC_Eigenvalues_Spanning_Tree_Number=-log(RUES2_ESC_Number_Non_Zero_Eigenvalues+1)+nansum(log(RUES2_ESC_Eigenvalues),2);
    RUES2_ESC_Eigenvalues_Spanning_Tree_Number(isnan(RUES2_ESC_Eigenvalues_Spanning_Tree_Number)) = 0;
    RUES2_ESC_Persistent_Spanning_Tree_Number=[RUES2_ESC_Persistent_Spanning_Tree_Number,RUES2_ESC_Eigenvalues_Spanning_Tree_Number];
end
RUES2_ESC_Average_Persistent_Multiplicity=mean(RUES2_ESC_Persistent_Multiplicity,2);
Persistent_Multiplicity=[Persistent_Multiplicity,RUES2_ESC_Average_Persistent_Multiplicity];

RUES2_ESC_Average_Persistent_Mean=mean(RUES2_ESC_Persistent_Mean,2);
Persistent_Mean=[Persistent_Mean,RUES2_ESC_Average_Persistent_Mean];

RUES2_ESC_Average_Persistent_Standard_Deviation=mean(RUES2_ESC_Persistent_Standard_Deviation,2);
Persistent_Standard_Deviation=[Persistent_Standard_Deviation,RUES2_ESC_Average_Persistent_Standard_Deviation];

RUES2_ESC_Average_Persistent_Maximum=mean(RUES2_ESC_Persistent_Maximum,2);
Persistent_Maximum=[Persistent_Maximum,RUES2_ESC_Average_Persistent_Maximum];

RUES2_ESC_Average_Persistent_Minimum=mean(RUES2_ESC_Persistent_Minimum,2);
Persistent_Minimum=[Persistent_Minimum,RUES2_ESC_Average_Persistent_Minimum];

RUES2_ESC_Average_Persistent_Laplacian_Graph_Energy=mean(RUES2_ESC_Persistent_Laplacian_Graph_Energy,2);
Persistent_Laplacian_Graph_Energy=[Persistent_Laplacian_Graph_Energy,RUES2_ESC_Average_Persistent_Laplacian_Graph_Energy];

RUES2_ESC_Average_Persistent_Generalized_Mean_Graph_Energy=mean(RUES2_ESC_Persistent_Generalized_Mean_Graph_Energy,2);
Persistent_Generalized_Mean_Graph_Energy=[Persistent_Generalized_Mean_Graph_Energy,RUES2_ESC_Average_Persistent_Generalized_Mean_Graph_Energy];

RUES2_ESC_Average_Persistent_Moment_Second_Order=mean(RUES2_ESC_Persistent_Moment_Second_Order,2);
Persistent_Moment_Second_Order=[Persistent_Moment_Second_Order,RUES2_ESC_Average_Persistent_Moment_Second_Order];

RUES2_ESC_Average_Persistent_Number_Non_Zero_Eigenvalue=mean(RUES2_ESC_Persistent_Number_Non_Zero_Eigenvalue,2);
Persistent_Number_Non_Zero_Eigenvalue=[Persistent_Number_Non_Zero_Eigenvalue,RUES2_ESC_Average_Persistent_Number_Non_Zero_Eigenvalue];

RUES2_ESC_Average_Persistent_Quasi_Wiener_Index=mean(RUES2_ESC_Persistent_Quasi_Wiener_Index,2);
Persistent_Quasi_Wiener_Index=[Persistent_Quasi_Wiener_Index,RUES2_ESC_Average_Persistent_Quasi_Wiener_Index];

RUES2_ESC_Average_Persistent_Spanning_Tree_Number=mean(RUES2_ESC_Persistent_Spanning_Tree_Number,2);
Persistent_Spanning_Tree_Number=[Persistent_Spanning_Tree_Number,RUES2_ESC_Average_Persistent_Spanning_Tree_Number];

% RUES2 FetalHeart
RUES2_FH_Persistent_Multiplicity=[];
RUES2_FH_Persistent_Mean=[];
RUES2_FH_Persistent_Standard_Deviation=[];
RUES2_FH_Persistent_Maximum=[];
RUES2_FH_Persistent_Minimum=[];
RUES2_FH_Persistent_Laplacian_Graph_Energy=[];
RUES2_FH_Persistent_Generalized_Mean_Graph_Energy=[];
RUES2_FH_Persistent_Moment_Second_Order=[];
RUES2_FH_Persistent_Number_Non_Zero_Eigenvalue=[];
RUES2_FH_Persistent_Quasi_Wiener_Index=[];
RUES2_FH_Persistent_Spanning_Tree_Number=[];
load('RUES2_FH_Chromosome_VR_L0_EV.mat');
RUES2_FH_Chr_Num=length(RUES2_FH_Chromosome_VR_L0_EV);
for RUES2_FH_Chr=1:RUES2_FH_Chr_Num
    RUES2_FH_Eigenvalues=RUES2_FH_Chromosome_VR_L0_EV{RUES2_FH_Chr};
    RUES2_FH_Eigenvalues(:,1)=[];
    RUES2_FH_Eigenvalues(101,:)=[];%Remove the eigenvalues when fully connected
    
    RUES2_FH_Number_Zero_Eigenvalues=sum(RUES2_FH_Eigenvalues==0,2);
    RUES2_FH_Persistent_Multiplicity=[RUES2_FH_Persistent_Multiplicity, RUES2_FH_Number_Zero_Eigenvalues];
    
    RUES2_FH_Eigenvalues(RUES2_FH_Eigenvalues==0) = NaN;
    
    RUES2_FH_Eigenvalues_Mean=nanmean(RUES2_FH_Eigenvalues,2);
    RUES2_FH_Eigenvalues_Mean(isnan(RUES2_FH_Eigenvalues_Mean)) = 0;
    RUES2_FH_Persistent_Mean=[RUES2_FH_Persistent_Mean,RUES2_FH_Eigenvalues_Mean];
    
    RUES2_FH_Eigenvalues_Standard_Deviation=nanstd(RUES2_FH_Eigenvalues,0,2);
    RUES2_FH_Eigenvalues_Standard_Deviation(isnan(RUES2_FH_Eigenvalues_Standard_Deviation)) = 0;
    RUES2_FH_Persistent_Standard_Deviation=[RUES2_FH_Persistent_Standard_Deviation,RUES2_FH_Eigenvalues_Standard_Deviation];
    
    RUES2_FH_Eigenvalues_Maximum=nanmax(RUES2_FH_Eigenvalues,[],2);
    RUES2_FH_Eigenvalues_Maximum(isnan(RUES2_FH_Eigenvalues_Maximum)) = 0;
    RUES2_FH_Persistent_Maximum=[RUES2_FH_Persistent_Maximum,RUES2_FH_Eigenvalues_Maximum];
    
    RUES2_FH_Eigenvalues_Minimum=nanmin(RUES2_FH_Eigenvalues,[],2);
    RUES2_FH_Eigenvalues_Minimum(isnan(RUES2_FH_Eigenvalues_Minimum)) = 0;
    RUES2_FH_Persistent_Minimum=[RUES2_FH_Persistent_Minimum,RUES2_FH_Eigenvalues_Minimum];
    
    RUES2_FH_Eigenvalues_Sum=nansum(RUES2_FH_Eigenvalues,2);
    RUES2_FH_Eigenvalues_Sum(isnan(RUES2_FH_Eigenvalues_Sum)) = 0;
    RUES2_FH_Persistent_Laplacian_Graph_Energy=[RUES2_FH_Persistent_Laplacian_Graph_Energy,RUES2_FH_Eigenvalues_Sum];
    
    RUES2_FH_Eigenvalues_Sum_Absolute_Deviation=nansum(abs(RUES2_FH_Eigenvalues-nanmean(RUES2_FH_Eigenvalues,2)),2);
    RUES2_FH_Eigenvalues_Sum_Absolute_Deviation(isnan(RUES2_FH_Eigenvalues_Sum_Absolute_Deviation)) = 0;
    RUES2_FH_Persistent_Generalized_Mean_Graph_Energy=[RUES2_FH_Persistent_Generalized_Mean_Graph_Energy,RUES2_FH_Eigenvalues_Sum_Absolute_Deviation];
    
    RUES2_FH_Eigenvalues_Moment_Second_Order=nansum(RUES2_FH_Eigenvalues.^2,2);
    RUES2_FH_Eigenvalues_Moment_Second_Order(isnan(RUES2_FH_Eigenvalues_Moment_Second_Order)) = 0;
    RUES2_FH_Persistent_Moment_Second_Order=[RUES2_FH_Persistent_Moment_Second_Order,RUES2_FH_Eigenvalues_Moment_Second_Order];
    
    RUES2_FH_Number_Non_Zero_Eigenvalues=sum(isnan(RUES2_FH_Eigenvalues)==0,2);
    RUES2_FH_Persistent_Number_Non_Zero_Eigenvalue=[RUES2_FH_Persistent_Number_Non_Zero_Eigenvalue,RUES2_FH_Number_Non_Zero_Eigenvalues];
    
    RUES2_FH_Eigenvalues_Quasi_Wiener_Index=nansum((RUES2_FH_Number_Non_Zero_Eigenvalues+1)./RUES2_FH_Eigenvalues,2);
    RUES2_FH_Eigenvalues_Quasi_Wiener_Index(isnan(RUES2_FH_Eigenvalues_Quasi_Wiener_Index)) = 0;
    RUES2_FH_Persistent_Quasi_Wiener_Index=[RUES2_FH_Persistent_Quasi_Wiener_Index,RUES2_FH_Eigenvalues_Quasi_Wiener_Index];
    
    RUES2_FH_Eigenvalues_Spanning_Tree_Number=-log(RUES2_FH_Number_Non_Zero_Eigenvalues+1)+nansum(log(RUES2_FH_Eigenvalues),2);
    RUES2_FH_Eigenvalues_Spanning_Tree_Number(isnan(RUES2_FH_Eigenvalues_Spanning_Tree_Number)) = 0;
    RUES2_FH_Persistent_Spanning_Tree_Number=[RUES2_FH_Persistent_Spanning_Tree_Number,RUES2_FH_Eigenvalues_Spanning_Tree_Number];
end
RUES2_FH_Average_Persistent_Multiplicity=mean(RUES2_FH_Persistent_Multiplicity,2);
Persistent_Multiplicity=[Persistent_Multiplicity,RUES2_FH_Average_Persistent_Multiplicity];

RUES2_FH_Average_Persistent_Mean=mean(RUES2_FH_Persistent_Mean,2);
Persistent_Mean=[Persistent_Mean,RUES2_FH_Average_Persistent_Mean];

RUES2_FH_Average_Persistent_Standard_Deviation=mean(RUES2_FH_Persistent_Standard_Deviation,2);
Persistent_Standard_Deviation=[Persistent_Standard_Deviation,RUES2_FH_Average_Persistent_Standard_Deviation];

RUES2_FH_Average_Persistent_Maximum=mean(RUES2_FH_Persistent_Maximum,2);
Persistent_Maximum=[Persistent_Maximum,RUES2_FH_Average_Persistent_Maximum];

RUES2_FH_Average_Persistent_Minimum=mean(RUES2_FH_Persistent_Minimum,2);
Persistent_Minimum=[Persistent_Minimum,RUES2_FH_Average_Persistent_Minimum];

RUES2_FH_Average_Persistent_Laplacian_Graph_Energy=mean(RUES2_FH_Persistent_Laplacian_Graph_Energy,2);
Persistent_Laplacian_Graph_Energy=[Persistent_Laplacian_Graph_Energy,RUES2_FH_Average_Persistent_Laplacian_Graph_Energy];

RUES2_FH_Average_Persistent_Generalized_Mean_Graph_Energy=mean(RUES2_FH_Persistent_Generalized_Mean_Graph_Energy,2);
Persistent_Generalized_Mean_Graph_Energy=[Persistent_Generalized_Mean_Graph_Energy,RUES2_FH_Average_Persistent_Generalized_Mean_Graph_Energy];

RUES2_FH_Average_Persistent_Moment_Second_Order=mean(RUES2_FH_Persistent_Moment_Second_Order,2);
Persistent_Moment_Second_Order=[Persistent_Moment_Second_Order,RUES2_FH_Average_Persistent_Moment_Second_Order];

RUES2_FH_Average_Persistent_Number_Non_Zero_Eigenvalue=mean(RUES2_FH_Persistent_Number_Non_Zero_Eigenvalue,2);
Persistent_Number_Non_Zero_Eigenvalue=[Persistent_Number_Non_Zero_Eigenvalue,RUES2_FH_Average_Persistent_Number_Non_Zero_Eigenvalue];

RUES2_FH_Average_Persistent_Quasi_Wiener_Index=mean(RUES2_FH_Persistent_Quasi_Wiener_Index,2);
Persistent_Quasi_Wiener_Index=[Persistent_Quasi_Wiener_Index,RUES2_FH_Average_Persistent_Quasi_Wiener_Index];

RUES2_FH_Average_Persistent_Spanning_Tree_Number=mean(RUES2_FH_Persistent_Spanning_Tree_Number,2);
Persistent_Spanning_Tree_Number=[Persistent_Spanning_Tree_Number,RUES2_FH_Average_Persistent_Spanning_Tree_Number];

% RUES2 MES
RUES2_MES_Persistent_Multiplicity=[];
RUES2_MES_Persistent_Mean=[];
RUES2_MES_Persistent_Standard_Deviation=[];
RUES2_MES_Persistent_Maximum=[];
RUES2_MES_Persistent_Minimum=[];
RUES2_MES_Persistent_Laplacian_Graph_Energy=[];
RUES2_MES_Persistent_Generalized_Mean_Graph_Energy=[];
RUES2_MES_Persistent_Moment_Second_Order=[];
RUES2_MES_Persistent_Number_Non_Zero_Eigenvalue=[];
RUES2_MES_Persistent_Quasi_Wiener_Index=[];
RUES2_MES_Persistent_Spanning_Tree_Number=[];
load('RUES2_MES_Chromosome_VR_L0_EV.mat');
RUES2_MES_Chr_Num=length(RUES2_MES_Chromosome_VR_L0_EV);
for RUES2_MES_Chr=1:RUES2_MES_Chr_Num
    RUES2_MES_Eigenvalues=RUES2_MES_Chromosome_VR_L0_EV{RUES2_MES_Chr};
    RUES2_MES_Eigenvalues(:,1)=[];
    RUES2_MES_Eigenvalues(101,:)=[];%Remove the eigenvalues when fully connected
    
    RUES2_MES_Number_Zero_Eigenvalues=sum(RUES2_MES_Eigenvalues==0,2);
    RUES2_MES_Persistent_Multiplicity=[RUES2_MES_Persistent_Multiplicity, RUES2_MES_Number_Zero_Eigenvalues];
    
    RUES2_MES_Eigenvalues(RUES2_MES_Eigenvalues==0) = NaN;
    
    RUES2_MES_Eigenvalues_Mean=nanmean(RUES2_MES_Eigenvalues,2);
    RUES2_MES_Eigenvalues_Mean(isnan(RUES2_MES_Eigenvalues_Mean)) = 0;
    RUES2_MES_Persistent_Mean=[RUES2_MES_Persistent_Mean,RUES2_MES_Eigenvalues_Mean];
    
    RUES2_MES_Eigenvalues_Standard_Deviation=nanstd(RUES2_MES_Eigenvalues,0,2);
    RUES2_MES_Eigenvalues_Standard_Deviation(isnan(RUES2_MES_Eigenvalues_Standard_Deviation)) = 0;
    RUES2_MES_Persistent_Standard_Deviation=[RUES2_MES_Persistent_Standard_Deviation,RUES2_MES_Eigenvalues_Standard_Deviation];
    
    RUES2_MES_Eigenvalues_Maximum=nanmax(RUES2_MES_Eigenvalues,[],2);
    RUES2_MES_Eigenvalues_Maximum(isnan(RUES2_MES_Eigenvalues_Maximum)) = 0;
    RUES2_MES_Persistent_Maximum=[RUES2_MES_Persistent_Maximum,RUES2_MES_Eigenvalues_Maximum];
    
    RUES2_MES_Eigenvalues_Minimum=nanmin(RUES2_MES_Eigenvalues,[],2);
    RUES2_MES_Eigenvalues_Minimum(isnan(RUES2_MES_Eigenvalues_Minimum)) = 0;
    RUES2_MES_Persistent_Minimum=[RUES2_MES_Persistent_Minimum,RUES2_MES_Eigenvalues_Minimum];
    
    RUES2_MES_Eigenvalues_Sum=nansum(RUES2_MES_Eigenvalues,2);
    RUES2_MES_Eigenvalues_Sum(isnan(RUES2_MES_Eigenvalues_Sum)) = 0;
    RUES2_MES_Persistent_Laplacian_Graph_Energy=[RUES2_MES_Persistent_Laplacian_Graph_Energy,RUES2_MES_Eigenvalues_Sum];
    
    RUES2_MES_Eigenvalues_Sum_Absolute_Deviation=nansum(abs(RUES2_MES_Eigenvalues-nanmean(RUES2_MES_Eigenvalues,2)),2);
    RUES2_MES_Eigenvalues_Sum_Absolute_Deviation(isnan(RUES2_MES_Eigenvalues_Sum_Absolute_Deviation)) = 0;
    RUES2_MES_Persistent_Generalized_Mean_Graph_Energy=[RUES2_MES_Persistent_Generalized_Mean_Graph_Energy,RUES2_MES_Eigenvalues_Sum_Absolute_Deviation];
    
    RUES2_MES_Eigenvalues_Moment_Second_Order=nansum(RUES2_MES_Eigenvalues.^2,2);
    RUES2_MES_Eigenvalues_Moment_Second_Order(isnan(RUES2_MES_Eigenvalues_Moment_Second_Order)) = 0;
    RUES2_MES_Persistent_Moment_Second_Order=[RUES2_MES_Persistent_Moment_Second_Order,RUES2_MES_Eigenvalues_Moment_Second_Order];
    
    RUES2_MES_Number_Non_Zero_Eigenvalues=sum(isnan(RUES2_MES_Eigenvalues)==0,2);
    RUES2_MES_Persistent_Number_Non_Zero_Eigenvalue=[RUES2_MES_Persistent_Number_Non_Zero_Eigenvalue,RUES2_MES_Number_Non_Zero_Eigenvalues];
    
    RUES2_MES_Eigenvalues_Quasi_Wiener_Index=nansum((RUES2_MES_Number_Non_Zero_Eigenvalues+1)./RUES2_MES_Eigenvalues,2);
    RUES2_MES_Eigenvalues_Quasi_Wiener_Index(isnan(RUES2_MES_Eigenvalues_Quasi_Wiener_Index)) = 0;
    RUES2_MES_Persistent_Quasi_Wiener_Index=[RUES2_MES_Persistent_Quasi_Wiener_Index,RUES2_MES_Eigenvalues_Quasi_Wiener_Index];
    
    RUES2_MES_Eigenvalues_Spanning_Tree_Number=-log(RUES2_MES_Number_Non_Zero_Eigenvalues+1)+nansum(log(RUES2_MES_Eigenvalues),2);
    RUES2_MES_Eigenvalues_Spanning_Tree_Number(isnan(RUES2_MES_Eigenvalues_Spanning_Tree_Number)) = 0;
    RUES2_MES_Persistent_Spanning_Tree_Number=[RUES2_MES_Persistent_Spanning_Tree_Number,RUES2_MES_Eigenvalues_Spanning_Tree_Number];
end
RUES2_MES_Average_Persistent_Multiplicity=mean(RUES2_MES_Persistent_Multiplicity,2);
Persistent_Multiplicity=[Persistent_Multiplicity,RUES2_MES_Average_Persistent_Multiplicity];

RUES2_MES_Average_Persistent_Mean=mean(RUES2_MES_Persistent_Mean,2);
Persistent_Mean=[Persistent_Mean,RUES2_MES_Average_Persistent_Mean];

RUES2_MES_Average_Persistent_Standard_Deviation=mean(RUES2_MES_Persistent_Standard_Deviation,2);
Persistent_Standard_Deviation=[Persistent_Standard_Deviation,RUES2_MES_Average_Persistent_Standard_Deviation];

RUES2_MES_Average_Persistent_Maximum=mean(RUES2_MES_Persistent_Maximum,2);
Persistent_Maximum=[Persistent_Maximum,RUES2_MES_Average_Persistent_Maximum];

RUES2_MES_Average_Persistent_Minimum=mean(RUES2_MES_Persistent_Minimum,2);
Persistent_Minimum=[Persistent_Minimum,RUES2_MES_Average_Persistent_Minimum];

RUES2_MES_Average_Persistent_Laplacian_Graph_Energy=mean(RUES2_MES_Persistent_Laplacian_Graph_Energy,2);
Persistent_Laplacian_Graph_Energy=[Persistent_Laplacian_Graph_Energy,RUES2_MES_Average_Persistent_Laplacian_Graph_Energy];

RUES2_MES_Average_Persistent_Generalized_Mean_Graph_Energy=mean(RUES2_MES_Persistent_Generalized_Mean_Graph_Energy,2);
Persistent_Generalized_Mean_Graph_Energy=[Persistent_Generalized_Mean_Graph_Energy,RUES2_MES_Average_Persistent_Generalized_Mean_Graph_Energy];

RUES2_MES_Average_Persistent_Moment_Second_Order=mean(RUES2_MES_Persistent_Moment_Second_Order,2);
Persistent_Moment_Second_Order=[Persistent_Moment_Second_Order,RUES2_MES_Average_Persistent_Moment_Second_Order];

RUES2_MES_Average_Persistent_Number_Non_Zero_Eigenvalue=mean(RUES2_MES_Persistent_Number_Non_Zero_Eigenvalue,2);
Persistent_Number_Non_Zero_Eigenvalue=[Persistent_Number_Non_Zero_Eigenvalue,RUES2_MES_Average_Persistent_Number_Non_Zero_Eigenvalue];

RUES2_MES_Average_Persistent_Quasi_Wiener_Index=mean(RUES2_MES_Persistent_Quasi_Wiener_Index,2);
Persistent_Quasi_Wiener_Index=[Persistent_Quasi_Wiener_Index,RUES2_MES_Average_Persistent_Quasi_Wiener_Index];

RUES2_MES_Average_Persistent_Spanning_Tree_Number=mean(RUES2_MES_Persistent_Spanning_Tree_Number,2);
Persistent_Spanning_Tree_Number=[Persistent_Spanning_Tree_Number,RUES2_MES_Average_Persistent_Spanning_Tree_Number];

%% WTC11
% WTC11 CM
WTC11_CM_Persistent_Multiplicity=[];
WTC11_CM_Persistent_Mean=[];
WTC11_CM_Persistent_Standard_Deviation=[];
WTC11_CM_Persistent_Maximum=[];
WTC11_CM_Persistent_Minimum=[];
WTC11_CM_Persistent_Laplacian_Graph_Energy=[];
WTC11_CM_Persistent_Generalized_Mean_Graph_Energy=[];
WTC11_CM_Persistent_Moment_Second_Order=[];
WTC11_CM_Persistent_Number_Non_Zero_Eigenvalue=[];
WTC11_CM_Persistent_Quasi_Wiener_Index=[];
WTC11_CM_Persistent_Spanning_Tree_Number=[];
load('WTC11_CM_Chromosome_VR_L0_EV.mat');
WTC11_CM_Chr_Num=length(WTC11_CM_Chromosome_VR_L0_EV);
for WTC11_CM_Chr=1:WTC11_CM_Chr_Num
    WTC11_CM_Eigenvalues=WTC11_CM_Chromosome_VR_L0_EV{WTC11_CM_Chr};
    WTC11_CM_Eigenvalues(:,1)=[];
    WTC11_CM_Eigenvalues(101,:)=[];%Remove the eigenvalues when fully connected
    
    WTC11_CM_Number_Zero_Eigenvalues=sum(WTC11_CM_Eigenvalues==0,2);
    WTC11_CM_Persistent_Multiplicity=[WTC11_CM_Persistent_Multiplicity, WTC11_CM_Number_Zero_Eigenvalues];
    
    WTC11_CM_Eigenvalues(WTC11_CM_Eigenvalues==0) = NaN;
    
    WTC11_CM_Eigenvalues_Mean=nanmean(WTC11_CM_Eigenvalues,2);
    WTC11_CM_Eigenvalues_Mean(isnan(WTC11_CM_Eigenvalues_Mean)) = 0;
    WTC11_CM_Persistent_Mean=[WTC11_CM_Persistent_Mean,WTC11_CM_Eigenvalues_Mean];
    
    WTC11_CM_Eigenvalues_Standard_Deviation=nanstd(WTC11_CM_Eigenvalues,0,2);
    WTC11_CM_Eigenvalues_Standard_Deviation(isnan(WTC11_CM_Eigenvalues_Standard_Deviation)) = 0;
    WTC11_CM_Persistent_Standard_Deviation=[WTC11_CM_Persistent_Standard_Deviation,WTC11_CM_Eigenvalues_Standard_Deviation];
    
    WTC11_CM_Eigenvalues_Maximum=nanmax(WTC11_CM_Eigenvalues,[],2);
    WTC11_CM_Eigenvalues_Maximum(isnan(WTC11_CM_Eigenvalues_Maximum)) = 0;
    WTC11_CM_Persistent_Maximum=[WTC11_CM_Persistent_Maximum,WTC11_CM_Eigenvalues_Maximum];
    
    WTC11_CM_Eigenvalues_Minimum=nanmin(WTC11_CM_Eigenvalues,[],2);
    WTC11_CM_Eigenvalues_Minimum(isnan(WTC11_CM_Eigenvalues_Minimum)) = 0;
    WTC11_CM_Persistent_Minimum=[WTC11_CM_Persistent_Minimum,WTC11_CM_Eigenvalues_Minimum];
    
    WTC11_CM_Eigenvalues_Sum=nansum(WTC11_CM_Eigenvalues,2);
    WTC11_CM_Eigenvalues_Sum(isnan(WTC11_CM_Eigenvalues_Sum)) = 0;
    WTC11_CM_Persistent_Laplacian_Graph_Energy=[WTC11_CM_Persistent_Laplacian_Graph_Energy,WTC11_CM_Eigenvalues_Sum];
    
    WTC11_CM_Eigenvalues_Sum_Absolute_Deviation=nansum(abs(WTC11_CM_Eigenvalues-nanmean(WTC11_CM_Eigenvalues,2)),2);
    WTC11_CM_Eigenvalues_Sum_Absolute_Deviation(isnan(WTC11_CM_Eigenvalues_Sum_Absolute_Deviation)) = 0;
    WTC11_CM_Persistent_Generalized_Mean_Graph_Energy=[WTC11_CM_Persistent_Generalized_Mean_Graph_Energy,WTC11_CM_Eigenvalues_Sum_Absolute_Deviation];
    
    WTC11_CM_Eigenvalues_Moment_Second_Order=nansum(WTC11_CM_Eigenvalues.^2,2);
    WTC11_CM_Eigenvalues_Moment_Second_Order(isnan(WTC11_CM_Eigenvalues_Moment_Second_Order)) = 0;
    WTC11_CM_Persistent_Moment_Second_Order=[WTC11_CM_Persistent_Moment_Second_Order,WTC11_CM_Eigenvalues_Moment_Second_Order];
    
    WTC11_CM_Number_Non_Zero_Eigenvalues=sum(isnan(WTC11_CM_Eigenvalues)==0,2);
    WTC11_CM_Persistent_Number_Non_Zero_Eigenvalue=[WTC11_CM_Persistent_Number_Non_Zero_Eigenvalue,WTC11_CM_Number_Non_Zero_Eigenvalues];
    
    WTC11_CM_Eigenvalues_Quasi_Wiener_Index=nansum((WTC11_CM_Number_Non_Zero_Eigenvalues+1)./WTC11_CM_Eigenvalues,2);
    WTC11_CM_Eigenvalues_Quasi_Wiener_Index(isnan(WTC11_CM_Eigenvalues_Quasi_Wiener_Index)) = 0;
    WTC11_CM_Persistent_Quasi_Wiener_Index=[WTC11_CM_Persistent_Quasi_Wiener_Index,WTC11_CM_Eigenvalues_Quasi_Wiener_Index];
    
    WTC11_CM_Eigenvalues_Spanning_Tree_Number=-log(WTC11_CM_Number_Non_Zero_Eigenvalues+1)+nansum(log(WTC11_CM_Eigenvalues),2);
    WTC11_CM_Eigenvalues_Spanning_Tree_Number(isnan(WTC11_CM_Eigenvalues_Spanning_Tree_Number)) = 0;
    WTC11_CM_Persistent_Spanning_Tree_Number=[WTC11_CM_Persistent_Spanning_Tree_Number,WTC11_CM_Eigenvalues_Spanning_Tree_Number];
end
WTC11_CM_Average_Persistent_Multiplicity=mean(WTC11_CM_Persistent_Multiplicity,2);
Persistent_Multiplicity=[Persistent_Multiplicity,WTC11_CM_Average_Persistent_Multiplicity];

WTC11_CM_Average_Persistent_Mean=mean(WTC11_CM_Persistent_Mean,2);
Persistent_Mean=[Persistent_Mean,WTC11_CM_Average_Persistent_Mean];

WTC11_CM_Average_Persistent_Standard_Deviation=mean(WTC11_CM_Persistent_Standard_Deviation,2);
Persistent_Standard_Deviation=[Persistent_Standard_Deviation,WTC11_CM_Average_Persistent_Standard_Deviation];

WTC11_CM_Average_Persistent_Maximum=mean(WTC11_CM_Persistent_Maximum,2);
Persistent_Maximum=[Persistent_Maximum,WTC11_CM_Average_Persistent_Maximum];

WTC11_CM_Average_Persistent_Minimum=mean(WTC11_CM_Persistent_Minimum,2);
Persistent_Minimum=[Persistent_Minimum,WTC11_CM_Average_Persistent_Minimum];

WTC11_CM_Average_Persistent_Laplacian_Graph_Energy=mean(WTC11_CM_Persistent_Laplacian_Graph_Energy,2);
Persistent_Laplacian_Graph_Energy=[Persistent_Laplacian_Graph_Energy,WTC11_CM_Average_Persistent_Laplacian_Graph_Energy];

WTC11_CM_Average_Persistent_Generalized_Mean_Graph_Energy=mean(WTC11_CM_Persistent_Generalized_Mean_Graph_Energy,2);
Persistent_Generalized_Mean_Graph_Energy=[Persistent_Generalized_Mean_Graph_Energy,WTC11_CM_Average_Persistent_Generalized_Mean_Graph_Energy];

WTC11_CM_Average_Persistent_Moment_Second_Order=mean(WTC11_CM_Persistent_Moment_Second_Order,2);
Persistent_Moment_Second_Order=[Persistent_Moment_Second_Order,WTC11_CM_Average_Persistent_Moment_Second_Order];

WTC11_CM_Average_Persistent_Number_Non_Zero_Eigenvalue=mean(WTC11_CM_Persistent_Number_Non_Zero_Eigenvalue,2);
Persistent_Number_Non_Zero_Eigenvalue=[Persistent_Number_Non_Zero_Eigenvalue,WTC11_CM_Average_Persistent_Number_Non_Zero_Eigenvalue];

WTC11_CM_Average_Persistent_Quasi_Wiener_Index=mean(WTC11_CM_Persistent_Quasi_Wiener_Index,2);
Persistent_Quasi_Wiener_Index=[Persistent_Quasi_Wiener_Index,WTC11_CM_Average_Persistent_Quasi_Wiener_Index];

WTC11_CM_Average_Persistent_Spanning_Tree_Number=mean(WTC11_CM_Persistent_Spanning_Tree_Number,2);
Persistent_Spanning_Tree_Number=[Persistent_Spanning_Tree_Number,WTC11_CM_Average_Persistent_Spanning_Tree_Number];

% WTC11 CP
WTC11_CP_Persistent_Multiplicity=[];
WTC11_CP_Persistent_Mean=[];
WTC11_CP_Persistent_Standard_Deviation=[];
WTC11_CP_Persistent_Maximum=[];
WTC11_CP_Persistent_Minimum=[];
WTC11_CP_Persistent_Laplacian_Graph_Energy=[];
WTC11_CP_Persistent_Generalized_Mean_Graph_Energy=[];
WTC11_CP_Persistent_Moment_Second_Order=[];
WTC11_CP_Persistent_Number_Non_Zero_Eigenvalue=[];
WTC11_CP_Persistent_Quasi_Wiener_Index=[];
WTC11_CP_Persistent_Spanning_Tree_Number=[];
load('WTC11_CP_Chromosome_VR_L0_EV.mat');
WTC11_CP_Chr_Num=length(WTC11_CP_Chromosome_VR_L0_EV);
for WTC11_CP_Chr=1:WTC11_CP_Chr_Num
    WTC11_CP_Eigenvalues=WTC11_CP_Chromosome_VR_L0_EV{WTC11_CP_Chr};
    WTC11_CP_Eigenvalues(:,1)=[];
    WTC11_CP_Eigenvalues(101,:)=[];%Remove the eigenvalues when fully connected
    
    WTC11_CP_Number_Zero_Eigenvalues=sum(WTC11_CP_Eigenvalues==0,2);
    WTC11_CP_Persistent_Multiplicity=[WTC11_CP_Persistent_Multiplicity, WTC11_CP_Number_Zero_Eigenvalues];
    
    WTC11_CP_Eigenvalues(WTC11_CP_Eigenvalues==0) = NaN;
    
    WTC11_CP_Eigenvalues_Mean=nanmean(WTC11_CP_Eigenvalues,2);
    WTC11_CP_Eigenvalues_Mean(isnan(WTC11_CP_Eigenvalues_Mean)) = 0;
    WTC11_CP_Persistent_Mean=[WTC11_CP_Persistent_Mean,WTC11_CP_Eigenvalues_Mean];
    
    WTC11_CP_Eigenvalues_Standard_Deviation=nanstd(WTC11_CP_Eigenvalues,0,2);
    WTC11_CP_Eigenvalues_Standard_Deviation(isnan(WTC11_CP_Eigenvalues_Standard_Deviation)) = 0;
    WTC11_CP_Persistent_Standard_Deviation=[WTC11_CP_Persistent_Standard_Deviation,WTC11_CP_Eigenvalues_Standard_Deviation];
    
    WTC11_CP_Eigenvalues_Maximum=nanmax(WTC11_CP_Eigenvalues,[],2);
    WTC11_CP_Eigenvalues_Maximum(isnan(WTC11_CP_Eigenvalues_Maximum)) = 0;
    WTC11_CP_Persistent_Maximum=[WTC11_CP_Persistent_Maximum,WTC11_CP_Eigenvalues_Maximum];
    
    WTC11_CP_Eigenvalues_Minimum=nanmin(WTC11_CP_Eigenvalues,[],2);
    WTC11_CP_Eigenvalues_Minimum(isnan(WTC11_CP_Eigenvalues_Minimum)) = 0;
    WTC11_CP_Persistent_Minimum=[WTC11_CP_Persistent_Minimum,WTC11_CP_Eigenvalues_Minimum];
    
    WTC11_CP_Eigenvalues_Sum=nansum(WTC11_CP_Eigenvalues,2);
    WTC11_CP_Eigenvalues_Sum(isnan(WTC11_CP_Eigenvalues_Sum)) = 0;
    WTC11_CP_Persistent_Laplacian_Graph_Energy=[WTC11_CP_Persistent_Laplacian_Graph_Energy,WTC11_CP_Eigenvalues_Sum];
    
    WTC11_CP_Eigenvalues_Sum_Absolute_Deviation=nansum(abs(WTC11_CP_Eigenvalues-nanmean(WTC11_CP_Eigenvalues,2)),2);
    WTC11_CP_Eigenvalues_Sum_Absolute_Deviation(isnan(WTC11_CP_Eigenvalues_Sum_Absolute_Deviation)) = 0;
    WTC11_CP_Persistent_Generalized_Mean_Graph_Energy=[WTC11_CP_Persistent_Generalized_Mean_Graph_Energy,WTC11_CP_Eigenvalues_Sum_Absolute_Deviation];
    
    WTC11_CP_Eigenvalues_Moment_Second_Order=nansum(WTC11_CP_Eigenvalues.^2,2);
    WTC11_CP_Eigenvalues_Moment_Second_Order(isnan(WTC11_CP_Eigenvalues_Moment_Second_Order)) = 0;
    WTC11_CP_Persistent_Moment_Second_Order=[WTC11_CP_Persistent_Moment_Second_Order,WTC11_CP_Eigenvalues_Moment_Second_Order];
    
    WTC11_CP_Number_Non_Zero_Eigenvalues=sum(isnan(WTC11_CP_Eigenvalues)==0,2);
    WTC11_CP_Persistent_Number_Non_Zero_Eigenvalue=[WTC11_CP_Persistent_Number_Non_Zero_Eigenvalue,WTC11_CP_Number_Non_Zero_Eigenvalues];
    
    WTC11_CP_Eigenvalues_Quasi_Wiener_Index=nansum((WTC11_CP_Number_Non_Zero_Eigenvalues+1)./WTC11_CP_Eigenvalues,2);
    WTC11_CP_Eigenvalues_Quasi_Wiener_Index(isnan(WTC11_CP_Eigenvalues_Quasi_Wiener_Index)) = 0;
    WTC11_CP_Persistent_Quasi_Wiener_Index=[WTC11_CP_Persistent_Quasi_Wiener_Index,WTC11_CP_Eigenvalues_Quasi_Wiener_Index];
    
    WTC11_CP_Eigenvalues_Spanning_Tree_Number=-log(WTC11_CP_Number_Non_Zero_Eigenvalues+1)+nansum(log(WTC11_CP_Eigenvalues),2);
    WTC11_CP_Eigenvalues_Spanning_Tree_Number(isnan(WTC11_CP_Eigenvalues_Spanning_Tree_Number)) = 0;
    WTC11_CP_Persistent_Spanning_Tree_Number=[WTC11_CP_Persistent_Spanning_Tree_Number,WTC11_CP_Eigenvalues_Spanning_Tree_Number];
end
WTC11_CP_Average_Persistent_Multiplicity=mean(WTC11_CP_Persistent_Multiplicity,2);
Persistent_Multiplicity=[Persistent_Multiplicity,WTC11_CP_Average_Persistent_Multiplicity];

WTC11_CP_Average_Persistent_Mean=mean(WTC11_CP_Persistent_Mean,2);
Persistent_Mean=[Persistent_Mean,WTC11_CP_Average_Persistent_Mean];

WTC11_CP_Average_Persistent_Standard_Deviation=mean(WTC11_CP_Persistent_Standard_Deviation,2);
Persistent_Standard_Deviation=[Persistent_Standard_Deviation,WTC11_CP_Average_Persistent_Standard_Deviation];

WTC11_CP_Average_Persistent_Maximum=mean(WTC11_CP_Persistent_Maximum,2);
Persistent_Maximum=[Persistent_Maximum,WTC11_CP_Average_Persistent_Maximum];

WTC11_CP_Average_Persistent_Minimum=mean(WTC11_CP_Persistent_Minimum,2);
Persistent_Minimum=[Persistent_Minimum,WTC11_CP_Average_Persistent_Minimum];

WTC11_CP_Average_Persistent_Laplacian_Graph_Energy=mean(WTC11_CP_Persistent_Laplacian_Graph_Energy,2);
Persistent_Laplacian_Graph_Energy=[Persistent_Laplacian_Graph_Energy,WTC11_CP_Average_Persistent_Laplacian_Graph_Energy];

WTC11_CP_Average_Persistent_Generalized_Mean_Graph_Energy=mean(WTC11_CP_Persistent_Generalized_Mean_Graph_Energy,2);
Persistent_Generalized_Mean_Graph_Energy=[Persistent_Generalized_Mean_Graph_Energy,WTC11_CP_Average_Persistent_Generalized_Mean_Graph_Energy];

WTC11_CP_Average_Persistent_Moment_Second_Order=mean(WTC11_CP_Persistent_Moment_Second_Order,2);
Persistent_Moment_Second_Order=[Persistent_Moment_Second_Order,WTC11_CP_Average_Persistent_Moment_Second_Order];

WTC11_CP_Average_Persistent_Number_Non_Zero_Eigenvalue=mean(WTC11_CP_Persistent_Number_Non_Zero_Eigenvalue,2);
Persistent_Number_Non_Zero_Eigenvalue=[Persistent_Number_Non_Zero_Eigenvalue,WTC11_CP_Average_Persistent_Number_Non_Zero_Eigenvalue];

WTC11_CP_Average_Persistent_Quasi_Wiener_Index=mean(WTC11_CP_Persistent_Quasi_Wiener_Index,2);
Persistent_Quasi_Wiener_Index=[Persistent_Quasi_Wiener_Index,WTC11_CP_Average_Persistent_Quasi_Wiener_Index];

WTC11_CP_Average_Persistent_Spanning_Tree_Number=mean(WTC11_CP_Persistent_Spanning_Tree_Number,2);
Persistent_Spanning_Tree_Number=[Persistent_Spanning_Tree_Number,WTC11_CP_Average_Persistent_Spanning_Tree_Number];

% WTC11 MES
WTC11_MES_Persistent_Multiplicity=[];
WTC11_MES_Persistent_Mean=[];
WTC11_MES_Persistent_Standard_Deviation=[];
WTC11_MES_Persistent_Maximum=[];
WTC11_MES_Persistent_Minimum=[];
WTC11_MES_Persistent_Laplacian_Graph_Energy=[];
WTC11_MES_Persistent_Generalized_Mean_Graph_Energy=[];
WTC11_MES_Persistent_Moment_Second_Order=[];
WTC11_MES_Persistent_Number_Non_Zero_Eigenvalue=[];
WTC11_MES_Persistent_Quasi_Wiener_Index=[];
WTC11_MES_Persistent_Spanning_Tree_Number=[];
load('WTC11_MES_Chromosome_VR_L0_EV.mat');
WTC11_MES_Chr_Num=length(WTC11_MES_Chromosome_VR_L0_EV);
for WTC11_MES_Chr=1:WTC11_MES_Chr_Num
    WTC11_MES_Eigenvalues=WTC11_MES_Chromosome_VR_L0_EV{WTC11_MES_Chr};
    WTC11_MES_Eigenvalues(:,1)=[];
    WTC11_MES_Eigenvalues(101,:)=[];%Remove the eigenvalues when fully connected
    
    WTC11_MES_Number_Zero_Eigenvalues=sum(WTC11_MES_Eigenvalues==0,2);
    WTC11_MES_Persistent_Multiplicity=[WTC11_MES_Persistent_Multiplicity, WTC11_MES_Number_Zero_Eigenvalues];
    
    WTC11_MES_Eigenvalues(WTC11_MES_Eigenvalues==0) = NaN;
    
    WTC11_MES_Eigenvalues_Mean=nanmean(WTC11_MES_Eigenvalues,2);
    WTC11_MES_Eigenvalues_Mean(isnan(WTC11_MES_Eigenvalues_Mean)) = 0;
    WTC11_MES_Persistent_Mean=[WTC11_MES_Persistent_Mean,WTC11_MES_Eigenvalues_Mean];
    
    WTC11_MES_Eigenvalues_Standard_Deviation=nanstd(WTC11_MES_Eigenvalues,0,2);
    WTC11_MES_Eigenvalues_Standard_Deviation(isnan(WTC11_MES_Eigenvalues_Standard_Deviation)) = 0;
    WTC11_MES_Persistent_Standard_Deviation=[WTC11_MES_Persistent_Standard_Deviation,WTC11_MES_Eigenvalues_Standard_Deviation];
    
    WTC11_MES_Eigenvalues_Maximum=nanmax(WTC11_MES_Eigenvalues,[],2);
    WTC11_MES_Eigenvalues_Maximum(isnan(WTC11_MES_Eigenvalues_Maximum)) = 0;
    WTC11_MES_Persistent_Maximum=[WTC11_MES_Persistent_Maximum,WTC11_MES_Eigenvalues_Maximum];
    
    WTC11_MES_Eigenvalues_Minimum=nanmin(WTC11_MES_Eigenvalues,[],2);
    WTC11_MES_Eigenvalues_Minimum(isnan(WTC11_MES_Eigenvalues_Minimum)) = 0;
    WTC11_MES_Persistent_Minimum=[WTC11_MES_Persistent_Minimum,WTC11_MES_Eigenvalues_Minimum];
    
    WTC11_MES_Eigenvalues_Sum=nansum(WTC11_MES_Eigenvalues,2);
    WTC11_MES_Eigenvalues_Sum(isnan(WTC11_MES_Eigenvalues_Sum)) = 0;
    WTC11_MES_Persistent_Laplacian_Graph_Energy=[WTC11_MES_Persistent_Laplacian_Graph_Energy,WTC11_MES_Eigenvalues_Sum];
    
    WTC11_MES_Eigenvalues_Sum_Absolute_Deviation=nansum(abs(WTC11_MES_Eigenvalues-nanmean(WTC11_MES_Eigenvalues,2)),2);
    WTC11_MES_Eigenvalues_Sum_Absolute_Deviation(isnan(WTC11_MES_Eigenvalues_Sum_Absolute_Deviation)) = 0;
    WTC11_MES_Persistent_Generalized_Mean_Graph_Energy=[WTC11_MES_Persistent_Generalized_Mean_Graph_Energy,WTC11_MES_Eigenvalues_Sum_Absolute_Deviation];
    
    WTC11_MES_Eigenvalues_Moment_Second_Order=nansum(WTC11_MES_Eigenvalues.^2,2);
    WTC11_MES_Eigenvalues_Moment_Second_Order(isnan(WTC11_MES_Eigenvalues_Moment_Second_Order)) = 0;
    WTC11_MES_Persistent_Moment_Second_Order=[WTC11_MES_Persistent_Moment_Second_Order,WTC11_MES_Eigenvalues_Moment_Second_Order];
    
    WTC11_MES_Number_Non_Zero_Eigenvalues=sum(isnan(WTC11_MES_Eigenvalues)==0,2);
    WTC11_MES_Persistent_Number_Non_Zero_Eigenvalue=[WTC11_MES_Persistent_Number_Non_Zero_Eigenvalue,WTC11_MES_Number_Non_Zero_Eigenvalues];
    
    WTC11_MES_Eigenvalues_Quasi_Wiener_Index=nansum((WTC11_MES_Number_Non_Zero_Eigenvalues+1)./WTC11_MES_Eigenvalues,2);
    WTC11_MES_Eigenvalues_Quasi_Wiener_Index(isnan(WTC11_MES_Eigenvalues_Quasi_Wiener_Index)) = 0;
    WTC11_MES_Persistent_Quasi_Wiener_Index=[WTC11_MES_Persistent_Quasi_Wiener_Index,WTC11_MES_Eigenvalues_Quasi_Wiener_Index];
    
    WTC11_MES_Eigenvalues_Spanning_Tree_Number=-log(WTC11_MES_Number_Non_Zero_Eigenvalues+1)+nansum(log(WTC11_MES_Eigenvalues),2);
    WTC11_MES_Eigenvalues_Spanning_Tree_Number(isnan(WTC11_MES_Eigenvalues_Spanning_Tree_Number)) = 0;
    WTC11_MES_Persistent_Spanning_Tree_Number=[WTC11_MES_Persistent_Spanning_Tree_Number,WTC11_MES_Eigenvalues_Spanning_Tree_Number];
end
WTC11_MES_Average_Persistent_Multiplicity=mean(WTC11_MES_Persistent_Multiplicity,2);
Persistent_Multiplicity=[Persistent_Multiplicity,WTC11_MES_Average_Persistent_Multiplicity];

WTC11_MES_Average_Persistent_Mean=mean(WTC11_MES_Persistent_Mean,2);
Persistent_Mean=[Persistent_Mean,WTC11_MES_Average_Persistent_Mean];

WTC11_MES_Average_Persistent_Standard_Deviation=mean(WTC11_MES_Persistent_Standard_Deviation,2);
Persistent_Standard_Deviation=[Persistent_Standard_Deviation,WTC11_MES_Average_Persistent_Standard_Deviation];

WTC11_MES_Average_Persistent_Maximum=mean(WTC11_MES_Persistent_Maximum,2);
Persistent_Maximum=[Persistent_Maximum,WTC11_MES_Average_Persistent_Maximum];

WTC11_MES_Average_Persistent_Minimum=mean(WTC11_MES_Persistent_Minimum,2);
Persistent_Minimum=[Persistent_Minimum,WTC11_MES_Average_Persistent_Minimum];

WTC11_MES_Average_Persistent_Laplacian_Graph_Energy=mean(WTC11_MES_Persistent_Laplacian_Graph_Energy,2);
Persistent_Laplacian_Graph_Energy=[Persistent_Laplacian_Graph_Energy,WTC11_MES_Average_Persistent_Laplacian_Graph_Energy];

WTC11_MES_Average_Persistent_Generalized_Mean_Graph_Energy=mean(WTC11_MES_Persistent_Generalized_Mean_Graph_Energy,2);
Persistent_Generalized_Mean_Graph_Energy=[Persistent_Generalized_Mean_Graph_Energy,WTC11_MES_Average_Persistent_Generalized_Mean_Graph_Energy];

WTC11_MES_Average_Persistent_Moment_Second_Order=mean(WTC11_MES_Persistent_Moment_Second_Order,2);
Persistent_Moment_Second_Order=[Persistent_Moment_Second_Order,WTC11_MES_Average_Persistent_Moment_Second_Order];

WTC11_MES_Average_Persistent_Number_Non_Zero_Eigenvalue=mean(WTC11_MES_Persistent_Number_Non_Zero_Eigenvalue,2);
Persistent_Number_Non_Zero_Eigenvalue=[Persistent_Number_Non_Zero_Eigenvalue,WTC11_MES_Average_Persistent_Number_Non_Zero_Eigenvalue];

WTC11_MES_Average_Persistent_Quasi_Wiener_Index=mean(WTC11_MES_Persistent_Quasi_Wiener_Index,2);
Persistent_Quasi_Wiener_Index=[Persistent_Quasi_Wiener_Index,WTC11_MES_Average_Persistent_Quasi_Wiener_Index];

WTC11_MES_Average_Persistent_Spanning_Tree_Number=mean(WTC11_MES_Persistent_Spanning_Tree_Number,2);
Persistent_Spanning_Tree_Number=[Persistent_Spanning_Tree_Number,WTC11_MES_Average_Persistent_Spanning_Tree_Number];

% WTC11 PSC
WTC11_PSC_Persistent_Multiplicity=[];
WTC11_PSC_Persistent_Mean=[];
WTC11_PSC_Persistent_Standard_Deviation=[];
WTC11_PSC_Persistent_Maximum=[];
WTC11_PSC_Persistent_Minimum=[];
WTC11_PSC_Persistent_Laplacian_Graph_Energy=[];
WTC11_PSC_Persistent_Generalized_Mean_Graph_Energy=[];
WTC11_PSC_Persistent_Moment_Second_Order=[];
WTC11_PSC_Persistent_Number_Non_Zero_Eigenvalue=[];
WTC11_PSC_Persistent_Quasi_Wiener_Index=[];
WTC11_PSC_Persistent_Spanning_Tree_Number=[];
load('WTC11_PSC_Chromosome_VR_L0_EV.mat');
WTC11_PSC_Chr_Num=length(WTC11_PSC_Chromosome_VR_L0_EV);
for WTC11_PSC_Chr=1:WTC11_PSC_Chr_Num
    WTC11_PSC_Eigenvalues=WTC11_PSC_Chromosome_VR_L0_EV{WTC11_PSC_Chr};
    WTC11_PSC_Eigenvalues(:,1)=[];
    WTC11_PSC_Eigenvalues(101,:)=[];%Remove the eigenvalues when fully connected
    
    WTC11_PSC_Number_Zero_Eigenvalues=sum(WTC11_PSC_Eigenvalues==0,2);
    WTC11_PSC_Persistent_Multiplicity=[WTC11_PSC_Persistent_Multiplicity, WTC11_PSC_Number_Zero_Eigenvalues];
    
    WTC11_PSC_Eigenvalues(WTC11_PSC_Eigenvalues==0) = NaN;
    
    WTC11_PSC_Eigenvalues_Mean=nanmean(WTC11_PSC_Eigenvalues,2);
    WTC11_PSC_Eigenvalues_Mean(isnan(WTC11_PSC_Eigenvalues_Mean)) = 0;
    WTC11_PSC_Persistent_Mean=[WTC11_PSC_Persistent_Mean,WTC11_PSC_Eigenvalues_Mean];
    
    WTC11_PSC_Eigenvalues_Standard_Deviation=nanstd(WTC11_PSC_Eigenvalues,0,2);
    WTC11_PSC_Eigenvalues_Standard_Deviation(isnan(WTC11_PSC_Eigenvalues_Standard_Deviation)) = 0;
    WTC11_PSC_Persistent_Standard_Deviation=[WTC11_PSC_Persistent_Standard_Deviation,WTC11_PSC_Eigenvalues_Standard_Deviation];
    
    WTC11_PSC_Eigenvalues_Maximum=nanmax(WTC11_PSC_Eigenvalues,[],2);
    WTC11_PSC_Eigenvalues_Maximum(isnan(WTC11_PSC_Eigenvalues_Maximum)) = 0;
    WTC11_PSC_Persistent_Maximum=[WTC11_PSC_Persistent_Maximum,WTC11_PSC_Eigenvalues_Maximum];
    
    WTC11_PSC_Eigenvalues_Minimum=nanmin(WTC11_PSC_Eigenvalues,[],2);
    WTC11_PSC_Eigenvalues_Minimum(isnan(WTC11_PSC_Eigenvalues_Minimum)) = 0;
    WTC11_PSC_Persistent_Minimum=[WTC11_PSC_Persistent_Minimum,WTC11_PSC_Eigenvalues_Minimum];
    
    WTC11_PSC_Eigenvalues_Sum=nansum(WTC11_PSC_Eigenvalues,2);
    WTC11_PSC_Eigenvalues_Sum(isnan(WTC11_PSC_Eigenvalues_Sum)) = 0;
    WTC11_PSC_Persistent_Laplacian_Graph_Energy=[WTC11_PSC_Persistent_Laplacian_Graph_Energy,WTC11_PSC_Eigenvalues_Sum];
    
    WTC11_PSC_Eigenvalues_Sum_Absolute_Deviation=nansum(abs(WTC11_PSC_Eigenvalues-nanmean(WTC11_PSC_Eigenvalues,2)),2);
    WTC11_PSC_Eigenvalues_Sum_Absolute_Deviation(isnan(WTC11_PSC_Eigenvalues_Sum_Absolute_Deviation)) = 0;
    WTC11_PSC_Persistent_Generalized_Mean_Graph_Energy=[WTC11_PSC_Persistent_Generalized_Mean_Graph_Energy,WTC11_PSC_Eigenvalues_Sum_Absolute_Deviation];
    
    WTC11_PSC_Eigenvalues_Moment_Second_Order=nansum(WTC11_PSC_Eigenvalues.^2,2);
    WTC11_PSC_Eigenvalues_Moment_Second_Order(isnan(WTC11_PSC_Eigenvalues_Moment_Second_Order)) = 0;
    WTC11_PSC_Persistent_Moment_Second_Order=[WTC11_PSC_Persistent_Moment_Second_Order,WTC11_PSC_Eigenvalues_Moment_Second_Order];
    
    WTC11_PSC_Number_Non_Zero_Eigenvalues=sum(isnan(WTC11_PSC_Eigenvalues)==0,2);
    WTC11_PSC_Persistent_Number_Non_Zero_Eigenvalue=[WTC11_PSC_Persistent_Number_Non_Zero_Eigenvalue,WTC11_PSC_Number_Non_Zero_Eigenvalues];
    
    WTC11_PSC_Eigenvalues_Quasi_Wiener_Index=nansum((WTC11_PSC_Number_Non_Zero_Eigenvalues+1)./WTC11_PSC_Eigenvalues,2);
    WTC11_PSC_Eigenvalues_Quasi_Wiener_Index(isnan(WTC11_PSC_Eigenvalues_Quasi_Wiener_Index)) = 0;
    WTC11_PSC_Persistent_Quasi_Wiener_Index=[WTC11_PSC_Persistent_Quasi_Wiener_Index,WTC11_PSC_Eigenvalues_Quasi_Wiener_Index];
    
    WTC11_PSC_Eigenvalues_Spanning_Tree_Number=-log(WTC11_PSC_Number_Non_Zero_Eigenvalues+1)+nansum(log(WTC11_PSC_Eigenvalues),2);
    WTC11_PSC_Eigenvalues_Spanning_Tree_Number(isnan(WTC11_PSC_Eigenvalues_Spanning_Tree_Number)) = 0;
    WTC11_PSC_Persistent_Spanning_Tree_Number=[WTC11_PSC_Persistent_Spanning_Tree_Number,WTC11_PSC_Eigenvalues_Spanning_Tree_Number];
end
WTC11_PSC_Average_Persistent_Multiplicity=mean(WTC11_PSC_Persistent_Multiplicity,2);
Persistent_Multiplicity=[Persistent_Multiplicity,WTC11_PSC_Average_Persistent_Multiplicity];

WTC11_PSC_Average_Persistent_Mean=mean(WTC11_PSC_Persistent_Mean,2);
Persistent_Mean=[Persistent_Mean,WTC11_PSC_Average_Persistent_Mean];

WTC11_PSC_Average_Persistent_Standard_Deviation=mean(WTC11_PSC_Persistent_Standard_Deviation,2);
Persistent_Standard_Deviation=[Persistent_Standard_Deviation,WTC11_PSC_Average_Persistent_Standard_Deviation];

WTC11_PSC_Average_Persistent_Maximum=mean(WTC11_PSC_Persistent_Maximum,2);
Persistent_Maximum=[Persistent_Maximum,WTC11_PSC_Average_Persistent_Maximum];

WTC11_PSC_Average_Persistent_Minimum=mean(WTC11_PSC_Persistent_Minimum,2);
Persistent_Minimum=[Persistent_Minimum,WTC11_PSC_Average_Persistent_Minimum];

WTC11_PSC_Average_Persistent_Laplacian_Graph_Energy=mean(WTC11_PSC_Persistent_Laplacian_Graph_Energy,2);
Persistent_Laplacian_Graph_Energy=[Persistent_Laplacian_Graph_Energy,WTC11_PSC_Average_Persistent_Laplacian_Graph_Energy];

WTC11_PSC_Average_Persistent_Generalized_Mean_Graph_Energy=mean(WTC11_PSC_Persistent_Generalized_Mean_Graph_Energy,2);
Persistent_Generalized_Mean_Graph_Energy=[Persistent_Generalized_Mean_Graph_Energy,WTC11_PSC_Average_Persistent_Generalized_Mean_Graph_Energy];

WTC11_PSC_Average_Persistent_Moment_Second_Order=mean(WTC11_PSC_Persistent_Moment_Second_Order,2);
Persistent_Moment_Second_Order=[Persistent_Moment_Second_Order,WTC11_PSC_Average_Persistent_Moment_Second_Order];

WTC11_PSC_Average_Persistent_Number_Non_Zero_Eigenvalue=mean(WTC11_PSC_Persistent_Number_Non_Zero_Eigenvalue,2);
Persistent_Number_Non_Zero_Eigenvalue=[Persistent_Number_Non_Zero_Eigenvalue,WTC11_PSC_Average_Persistent_Number_Non_Zero_Eigenvalue];

WTC11_PSC_Average_Persistent_Quasi_Wiener_Index=mean(WTC11_PSC_Persistent_Quasi_Wiener_Index,2);
Persistent_Quasi_Wiener_Index=[Persistent_Quasi_Wiener_Index,WTC11_PSC_Average_Persistent_Quasi_Wiener_Index];

WTC11_PSC_Average_Persistent_Spanning_Tree_Number=mean(WTC11_PSC_Persistent_Spanning_Tree_Number,2);
Persistent_Spanning_Tree_Number=[Persistent_Spanning_Tree_Number,WTC11_PSC_Average_Persistent_Spanning_Tree_Number];

%% Structural classification (SC): row: 14, column: 11*24*100=26400.
%% Local attributes (LA): row: 14, column: 100
for i=1:24
    
    SC_H1_ESC_PD=[];
    SC_H1_ME_PD=[];
    SC_H1_MS_PD=[];
    SC_H1_NP_PD=[];
    SC_H1_TB_PD=[];
    SC_RUES2_CM_PD=[];
    SC_RUES2_CP_PD=[];
    SC_RUES2_ESC_PD=[];
    SC_RUES2_FH_PD=[];
    SC_RUES2_MES_PD=[];
    SC_WTC11_CM_PD=[];
    SC_WTC11_CP_PD=[];
    SC_WTC11_MES_PD=[];
    SC_WTC11_PSC_PD=[];
    
    SC_Persistent_Descriptors=[];
    
    SC_H1_ESC_PD=[SC_H1_ESC_PD ...
         H1_ESC_Persistent_Generalized_Mean_Graph_Energy(:,i) ...
        % H1_ESC_Persistent_Laplacian_Graph_Energy(:,i) ...
        % H1_ESC_Persistent_Maximum(:,i) ...
        % H1_ESC_Persistent_Mean(:,i) ...
        % H1_ESC_Persistent_Minimum(:,i) ...
        % H1_ESC_Persistent_Moment_Second_Order(:,i) ...
        % H1_ESC_Persistent_Multiplicity(:,i) ...
        % H1_ESC_Persistent_Number_Non_Zero_Eigenvalue(:,i) ...
        % H1_ESC_Persistent_Quasi_Wiener_Index(:,i) ...
        % H1_ESC_Persistent_Spanning_Tree_Number(:,i) ...
        % H1_ESC_Persistent_Standard_Deviation(:,i)...
        ];
    SC_H1_ME_PD=[SC_H1_ME_PD ...
        H1_ME_Persistent_Generalized_Mean_Graph_Energy(:,i) ...
        % H1_ME_Persistent_Laplacian_Graph_Energy(:,i) ...
        % H1_ME_Persistent_Maximum(:,i) ...
        % H1_ME_Persistent_Mean(:,i) ...
        % H1_ME_Persistent_Minimum(:,i) ...
        % H1_ME_Persistent_Moment_Second_Order(:,i) ...
        % H1_ME_Persistent_Multiplicity(:,i) ...
        % H1_ME_Persistent_Number_Non_Zero_Eigenvalue(:,i) ...
        % H1_ME_Persistent_Quasi_Wiener_Index(:,i) ...
        % H1_ME_Persistent_Spanning_Tree_Number(:,i) ...
        % H1_ME_Persistent_Standard_Deviation(:,i) ...
        ];
    SC_H1_MS_PD=[SC_H1_MS_PD ...
        H1_MS_Persistent_Generalized_Mean_Graph_Energy(:,i) ...
        % H1_MS_Persistent_Laplacian_Graph_Energy(:,i) ...
        % H1_MS_Persistent_Maximum(:,i) ...
        % H1_MS_Persistent_Mean(:,i) ...
        % H1_MS_Persistent_Minimum(:,i) ...
        % H1_MS_Persistent_Moment_Second_Order(:,i) ...
        % H1_MS_Persistent_Multiplicity(:,i) ...
        % H1_MS_Persistent_Number_Non_Zero_Eigenvalue(:,i) ...
        % H1_MS_Persistent_Quasi_Wiener_Index(:,i) ...
        % H1_MS_Persistent_Spanning_Tree_Number(:,i) ...
        % H1_MS_Persistent_Standard_Deviation(:,i) ...
        ];
    SC_H1_NP_PD=[SC_H1_NP_PD ...
        H1_NP_Persistent_Generalized_Mean_Graph_Energy(:,i) ...
        % H1_NP_Persistent_Laplacian_Graph_Energy(:,i) ...
        % H1_NP_Persistent_Maximum(:,i) ...
        % H1_NP_Persistent_Mean(:,i) ...
        % H1_NP_Persistent_Minimum(:,i) ...
        % H1_NP_Persistent_Moment_Second_Order(:,i) ...
        % H1_NP_Persistent_Multiplicity(:,i) ...
        % H1_NP_Persistent_Number_Non_Zero_Eigenvalue(:,i) ...
        % H1_NP_Persistent_Quasi_Wiener_Index(:,i) ...
        % H1_NP_Persistent_Spanning_Tree_Number(:,i) ...
        % H1_NP_Persistent_Standard_Deviation(:,i) ...
        ];
    SC_H1_TB_PD=[SC_H1_TB_PD ...
        H1_TB_Persistent_Generalized_Mean_Graph_Energy(:,i) ...
        % H1_TB_Persistent_Laplacian_Graph_Energy(:,i) ...
        % H1_TB_Persistent_Maximum(:,i) ...
        % H1_TB_Persistent_Mean(:,i) ...
        % H1_TB_Persistent_Minimum(:,i) ...
        % H1_TB_Persistent_Moment_Second_Order(:,i) ...
        % H1_TB_Persistent_Multiplicity(:,i) ...
        % H1_TB_Persistent_Number_Non_Zero_Eigenvalue(:,i) ...
        % H1_TB_Persistent_Quasi_Wiener_Index(:,i) ...
        % H1_TB_Persistent_Spanning_Tree_Number(:,i) ...
        % H1_TB_Persistent_Standard_Deviation(:,i) ...
        ];
    SC_RUES2_CM_PD=[SC_RUES2_CM_PD ...
        RUES2_CM_Persistent_Generalized_Mean_Graph_Energy(:,i) ...
        % RUES2_CM_Persistent_Laplacian_Graph_Energy(:,i) ...
        % RUES2_CM_Persistent_Maximum(:,i) ...
        % RUES2_CM_Persistent_Mean(:,i) ...
        % RUES2_CM_Persistent_Minimum(:,i) ...
        % RUES2_CM_Persistent_Moment_Second_Order(:,i) ...
        % RUES2_CM_Persistent_Multiplicity(:,i) ...
        % RUES2_CM_Persistent_Number_Non_Zero_Eigenvalue(:,i) ...
        % RUES2_CM_Persistent_Quasi_Wiener_Index(:,i) ...
        % RUES2_CM_Persistent_Spanning_Tree_Number(:,i) ...
        % RUES2_CM_Persistent_Standard_Deviation(:,i) ...
        ];
    SC_RUES2_CP_PD=[SC_RUES2_CP_PD ...
        RUES2_CP_Persistent_Generalized_Mean_Graph_Energy(:,i) ...
        % RUES2_CP_Persistent_Laplacian_Graph_Energy(:,i) ...
        % RUES2_CP_Persistent_Maximum(:,i) ...
        % RUES2_CP_Persistent_Mean(:,i) ...
        % RUES2_CP_Persistent_Minimum(:,i) ...
        % RUES2_CP_Persistent_Moment_Second_Order(:,i) ...
        % RUES2_CP_Persistent_Multiplicity(:,i) ...
        % RUES2_CP_Persistent_Number_Non_Zero_Eigenvalue(:,i) ...
        % RUES2_CP_Persistent_Quasi_Wiener_Index(:,i) ...
        % RUES2_CP_Persistent_Spanning_Tree_Number(:,i) ...
        % RUES2_CP_Persistent_Standard_Deviation(:,i) ...
        ];
    SC_RUES2_ESC_PD=[SC_RUES2_ESC_PD ...
        RUES2_ESC_Persistent_Generalized_Mean_Graph_Energy(:,i) ...
        % RUES2_ESC_Persistent_Laplacian_Graph_Energy(:,i) ...
        % RUES2_ESC_Persistent_Maximum(:,i) ...
        % RUES2_ESC_Persistent_Mean(:,i) ...
        % RUES2_ESC_Persistent_Minimum(:,i) ...
        % RUES2_ESC_Persistent_Moment_Second_Order(:,i) ...
        % RUES2_ESC_Persistent_Multiplicity(:,i) ...
        % RUES2_ESC_Persistent_Number_Non_Zero_Eigenvalue(:,i) ...
        % RUES2_ESC_Persistent_Quasi_Wiener_Index(:,i) ...
        % RUES2_ESC_Persistent_Spanning_Tree_Number(:,i) ...
        % RUES2_ESC_Persistent_Standard_Deviation(:,i) ...
        ];
    SC_RUES2_FH_PD=[SC_RUES2_FH_PD ...
        RUES2_FH_Persistent_Generalized_Mean_Graph_Energy(:,i) ...
        % RUES2_FH_Persistent_Laplacian_Graph_Energy(:,i) ...
        % RUES2_FH_Persistent_Maximum(:,i) ...
        % RUES2_FH_Persistent_Mean(:,i) ...
        % RUES2_FH_Persistent_Minimum(:,i) ...
        % RUES2_FH_Persistent_Moment_Second_Order(:,i) ...
        % RUES2_FH_Persistent_Multiplicity(:,i) ...
        % RUES2_FH_Persistent_Number_Non_Zero_Eigenvalue(:,i) ...
        % RUES2_FH_Persistent_Quasi_Wiener_Index(:,i) ...
        % RUES2_FH_Persistent_Spanning_Tree_Number(:,i) ...
        % RUES2_FH_Persistent_Standard_Deviation(:,i) ...
        ];
    SC_RUES2_MES_PD=[SC_RUES2_MES_PD ...
        RUES2_MES_Persistent_Generalized_Mean_Graph_Energy(:,i) ...
        % RUES2_MES_Persistent_Laplacian_Graph_Energy(:,i) ...
        % RUES2_MES_Persistent_Maximum(:,i) ...
        % RUES2_MES_Persistent_Mean(:,i) ...
        % RUES2_MES_Persistent_Minimum(:,i) ...
        % RUES2_MES_Persistent_Moment_Second_Order(:,i) ...
        % RUES2_MES_Persistent_Multiplicity(:,i) ...
        % RUES2_MES_Persistent_Number_Non_Zero_Eigenvalue(:,i) ...
        % RUES2_MES_Persistent_Quasi_Wiener_Index(:,i) ...
        % RUES2_MES_Persistent_Spanning_Tree_Number(:,i) ...
        % RUES2_MES_Persistent_Standard_Deviation(:,i) ...
        ];
    SC_WTC11_CM_PD=[SC_WTC11_CM_PD ...
        WTC11_CM_Persistent_Generalized_Mean_Graph_Energy(:,i) ...
        % WTC11_CM_Persistent_Laplacian_Graph_Energy(:,i) ...
        % WTC11_CM_Persistent_Maximum(:,i) ...
        % WTC11_CM_Persistent_Mean(:,i) ...
        % WTC11_CM_Persistent_Minimum(:,i) ...
        % WTC11_CM_Persistent_Moment_Second_Order(:,i) ...
        % WTC11_CM_Persistent_Multiplicity(:,i) ...
        % WTC11_CM_Persistent_Number_Non_Zero_Eigenvalue(:,i) ...
        % WTC11_CM_Persistent_Quasi_Wiener_Index(:,i) ...
        % WTC11_CM_Persistent_Spanning_Tree_Number(:,i) ...
        % WTC11_CM_Persistent_Standard_Deviation(:,i) ...
        ];
    SC_WTC11_CP_PD=[SC_WTC11_CP_PD ...
        WTC11_CP_Persistent_Generalized_Mean_Graph_Energy(:,i) ...
        % WTC11_CP_Persistent_Laplacian_Graph_Energy(:,i) ...
        % WTC11_CP_Persistent_Maximum(:,i) ...
        % WTC11_CP_Persistent_Mean(:,i) ...
        % WTC11_CP_Persistent_Minimum(:,i) ...
        % WTC11_CP_Persistent_Moment_Second_Order(:,i) ...
        % WTC11_CP_Persistent_Multiplicity(:,i) ...
        % WTC11_CP_Persistent_Number_Non_Zero_Eigenvalue(:,i) ...
        % WTC11_CP_Persistent_Quasi_Wiener_Index(:,i) ...
        % WTC11_CP_Persistent_Spanning_Tree_Number(:,i) ...
        % WTC11_CP_Persistent_Standard_Deviation(:,i) ...
        ];
    SC_WTC11_MES_PD=[SC_WTC11_MES_PD ...
        WTC11_MES_Persistent_Generalized_Mean_Graph_Energy(:,i) ...
        % WTC11_MES_Persistent_Laplacian_Graph_Energy(:,i) ...
        % WTC11_MES_Persistent_Maximum(:,i) ...
        % WTC11_MES_Persistent_Mean(:,i) ...
        % WTC11_MES_Persistent_Minimum(:,i) ...
        % WTC11_MES_Persistent_Moment_Second_Order(:,i) ...
        % WTC11_MES_Persistent_Multiplicity(:,i) ...
        % WTC11_MES_Persistent_Number_Non_Zero_Eigenvalue(:,i) ...
        % WTC11_MES_Persistent_Quasi_Wiener_Index(:,i) ...
        % WTC11_MES_Persistent_Spanning_Tree_Number(:,i) ...
        % WTC11_MES_Persistent_Standard_Deviation(:,i) ...
        ];
    SC_WTC11_PSC_PD=[SC_WTC11_PSC_PD ...
        WTC11_PSC_Persistent_Generalized_Mean_Graph_Energy(:,i) ...
        % WTC11_PSC_Persistent_Laplacian_Graph_Energy(:,i) ...
        % WTC11_PSC_Persistent_Maximum(:,i) ...
        % WTC11_PSC_Persistent_Mean(:,i) ...
        % WTC11_PSC_Persistent_Minimum(:,i) ...
        % WTC11_PSC_Persistent_Moment_Second_Order(:,i) ...
        % WTC11_PSC_Persistent_Multiplicity(:,i) ...
        % WTC11_PSC_Persistent_Number_Non_Zero_Eigenvalue(:,i) ...
        % WTC11_PSC_Persistent_Quasi_Wiener_Index(:,i) ...
        % WTC11_PSC_Persistent_Spanning_Tree_Number(:,i) ...
        % WTC11_PSC_Persistent_Standard_Deviation(:,i) ...
        ];
    SC_Persistent_Descriptors=[SC_H1_ESC_PD';SC_H1_ME_PD';SC_H1_MS_PD';SC_H1_NP_PD';SC_H1_TB_PD';...
        SC_RUES2_CM_PD';SC_RUES2_CP_PD';SC_RUES2_ESC_PD';SC_RUES2_FH_PD';SC_RUES2_MES_PD';...
        SC_WTC11_CM_PD';SC_WTC11_CP_PD';SC_WTC11_MES_PD';SC_WTC11_PSC_PD'];
    
    %% t-SNE
    %t-SNE1 t-SNE2
    figure('Color','white','Position',[400 0 1000 1000]);
    subplot(2,1,1);
    [SC_Y_t_SNE_2,SC_loss_2] = tsne(SC_Persistent_Descriptors,'Algorithm','exact','Distance','cosine','Perplexity',6);
    SC_Persistent_Descriptors=SC_Persistent_Descriptors';
    SC_Persistent_Generalized_Mean_Graph_Energy_Attributes{i}=SC_Persistent_Descriptors;
    % SC_Persistent_Laplacian_Graph_Energy{i}=SC_Persistent_Descriptors;
    % SC_Persistent_Maximum{i}=SC_Persistent_Descriptors;
    % SC_Persistent_Mean{i}=SC_Persistent_Descriptors;
    % SC_Persistent_Minimum{i}=SC_Persistent_Descriptors;
    % SC_Persistent_Moment_Second_Order{i}=SC_Persistent_Descriptors;
    % SC_Persistent_Multiplicity{i}=SC_Persistent_Descriptors;
    % SC_Persistent_Number_Non_Zero_Eigenvalue{i}=SC_Persistent_Descriptors;
    % SC_Persistent_Quasi_Wiener_Index{i}=SC_Persistent_Descriptors;
    % SC_Persistent_Spanning_Tree_Number{i}=SC_Persistent_Descriptors;
    % SC_Persistent_Standard_Deviation{i}=SC_Persistent_Descriptors;
    scatter(SC_Y_t_SNE_2(1,1),SC_Y_t_SNE_2(1,2),'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','r','LineWidth',1);%H1 ESC
    text(SC_Y_t_SNE_2(1,1),SC_Y_t_SNE_2(1,2),'H1 ESC','FontSize',12,'FontName','Times New Roman','HorizontalAlignment','center','VerticalAlignment','bottom');
    hold on
    scatter(SC_Y_t_SNE_2(2,1),SC_Y_t_SNE_2(2,2),'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','r','LineWidth',1);%H1 ME
    text(SC_Y_t_SNE_2(2,1),SC_Y_t_SNE_2(2,2),'H1 ME','FontSize',12,'FontName','Times New Roman','HorizontalAlignment','center','VerticalAlignment','bottom');
    hold on
    scatter(SC_Y_t_SNE_2(3,1),SC_Y_t_SNE_2(3,2),'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','r','LineWidth',1);%H1 MS
    text(SC_Y_t_SNE_2(3,1),SC_Y_t_SNE_2(3,2),'H1 MS','FontSize',12,'FontName','Times New Roman','HorizontalAlignment','center','VerticalAlignment','bottom');
    hold on
    scatter(SC_Y_t_SNE_2(4,1),SC_Y_t_SNE_2(4,2),'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','r','LineWidth',1);%H1 NP
    text(SC_Y_t_SNE_2(4,1),SC_Y_t_SNE_2(4,2),'H1 NP','FontSize',12,'FontName','Times New Roman','HorizontalAlignment','center','VerticalAlignment','bottom');
    hold on
    scatter(SC_Y_t_SNE_2(5,1),SC_Y_t_SNE_2(5,2),'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','r','LineWidth',1);%H1 TB
    text(SC_Y_t_SNE_2(5,1),SC_Y_t_SNE_2(5,2),'H1 TB','FontSize',12,'FontName','Times New Roman','HorizontalAlignment','center','VerticalAlignment','bottom');
    hold on
    scatter(SC_Y_t_SNE_2(6,1),SC_Y_t_SNE_2(6,2),'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','b','LineWidth',1);%RUES2 CM
    text(SC_Y_t_SNE_2(6,1),SC_Y_t_SNE_2(6,2),'RUES2 CM','FontSize',12,'FontName','Times New Roman','HorizontalAlignment','center','VerticalAlignment','bottom');
    hold on
    scatter(SC_Y_t_SNE_2(7,1),SC_Y_t_SNE_2(7,2),'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','b','LineWidth',1);%RUES2 CP
    text(SC_Y_t_SNE_2(7,1),SC_Y_t_SNE_2(7,2),'RUES2 CP','FontSize',12,'FontName','Times New Roman','HorizontalAlignment','center','VerticalAlignment','bottom');
    hold on
    scatter(SC_Y_t_SNE_2(8,1),SC_Y_t_SNE_2(8,2),'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','b','LineWidth',1);%RUES2 ESC
    text(SC_Y_t_SNE_2(8,1),SC_Y_t_SNE_2(8,2),'RUES2 ESC','FontSize',12,'FontName','Times New Roman','HorizontalAlignment','center','VerticalAlignment','bottom');
    hold on
    scatter(SC_Y_t_SNE_2(9,1),SC_Y_t_SNE_2(9,2),'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','b','LineWidth',1);%RUES2 FetalHeart
    text(SC_Y_t_SNE_2(9,1),SC_Y_t_SNE_2(9,2),'RUES2 FetalHeart','FontSize',12,'FontName','Times New Roman','HorizontalAlignment','center','VerticalAlignment','bottom');
    hold on
    scatter(SC_Y_t_SNE_2(10,1),SC_Y_t_SNE_2(10,2),'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','b','LineWidth',1);%RUES2 MES
    text(SC_Y_t_SNE_2(10,1),SC_Y_t_SNE_2(10,2),'RUES2 MES','FontSize',12,'FontName','Times New Roman','HorizontalAlignment','center','VerticalAlignment','bottom');
    hold on
    scatter(SC_Y_t_SNE_2(11,1),SC_Y_t_SNE_2(11,2),'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','k','LineWidth',1);%WTC11 CM
    text(SC_Y_t_SNE_2(11,1),SC_Y_t_SNE_2(11,2),'WTC11 CM','FontSize',12,'FontName','Times New Roman','HorizontalAlignment','center','VerticalAlignment','bottom');
    hold on
    scatter(SC_Y_t_SNE_2(12,1),SC_Y_t_SNE_2(12,2),'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','k','LineWidth',1);%WTC11 CP
    text(SC_Y_t_SNE_2(12,1),SC_Y_t_SNE_2(12,2),'WTC11 CP','FontSize',12,'FontName','Times New Roman','HorizontalAlignment','center','VerticalAlignment','bottom');
    hold on
    scatter(SC_Y_t_SNE_2(13,1),SC_Y_t_SNE_2(13,2),'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','k','LineWidth',1);%WTC11 MES
    text(SC_Y_t_SNE_2(13,1),SC_Y_t_SNE_2(13,2),'WTC11 MES','FontSize',12,'FontName','Times New Roman','HorizontalAlignment','center','VerticalAlignment','bottom');
    hold on
    scatter(SC_Y_t_SNE_2(14,1),SC_Y_t_SNE_2(14,2),'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','k','LineWidth',1);%WTC11 PSC
    text(SC_Y_t_SNE_2(14,1),SC_Y_t_SNE_2(14,2),'WTC11 PSC','FontSize',12,'FontName','Times New Roman','HorizontalAlignment','center','VerticalAlignment','bottom');
    hold off
    title('t-SNE','FontSize',20,'FontName','Times New Roman');
    SC_Persistent_Generalized_Mean_Graph_Energy_loss{i}=SC_loss_2;
    % SC_Persistent_Laplacian_Graph_Energy_loss{i}=SC_loss_2;
    % SC_Persistent_Maximum_loss{i}=SC_loss_2;
    % SC_Persistent_Mean_loss{i}=SC_loss_2;
    % SC_Persistent_Minimum_loss{i}=SC_loss_2;
    % SC_Persistent_Moment_Second_Order_loss{i}=SC_loss_2;
    % SC_Persistent_Multiplicity_loss{i}=SC_loss_2;
    % SC_Persistent_Number_Non_Zero_Eigenvalue_loss{i}=SC_loss_2;
    % SC_Persistent_Quasi_Wiener_Index_loss{i}=SC_loss_2;
    % SC_Persistent_Spanning_Tree_Number_loss{i}=SC_loss_2;
    % SC_Persistent_Standard_Deviation_loss{i}=SC_loss_2;
    SC_Y_t_SNE=reshape(SC_Y_t_SNE_2',1,28);
    SC_Persistent_Generalized_Mean_Graph_Energy_Y_t_SNE{i}=SC_Y_t_SNE;
    % SC_Persistent_Laplacian_Graph_Energy_Y_t_SNE{i}=SC_Y_t_SNE;
    % SC_Persistent_Maximum_Y_t_SNE{i}=SC_Y_t_SNE;
    % SC_Persistent_Mean_Y_t_SNE{i}=SC_Y_t_SNE;
    % SC_Persistent_Minimum_Y_t_SNE{i}=SC_Y_t_SNE;
    % SC_Persistent_Moment_Second_Order_Y_t_SNE{i}=SC_Y_t_SNE;
    % SC_Persistent_Multiplicity_Y_t_SNE{i}=SC_Y_t_SNE;
    % SC_Persistent_Number_Non_Zero_Eigenvalue_Y_t_SNE{i}=SC_Y_t_SNE;
    % SC_Persistent_Quasi_Wiener_Index_Y_t_SNE{i}=SC_Y_t_SNE;
    % SC_Persistent_Spanning_Tree_Number_Y_t_SNE{i}=SC_Y_t_SNE;
    % SC_Persistent_Standard_Deviation_Y_t_SNE{i}=SC_Y_t_SNE;
    %k-means
    [SC_idx_2,SC_Centroids_2] = kmeans(SC_Y_t_SNE_2,3);
    SC_Persistent_Generalized_Mean_Graph_Energy_idx_2{i}=SC_idx_2';
    % SC_Persistent_Laplacian_Graph_Energy_idx_2{i}=SC_idx_2';
    % SC_Persistent_Maximum_idx_2{i}=SC_idx_2';
    % SC_Persistent_Mean_idx_2{i}=SC_idx_2';
    % SC_Persistent_Minimum_idx_2{i}=SC_idx_2';
    % SC_Persistent_Moment_Second_Order_idx_2{i}=SC_idx_2';
    % SC_Persistent_Multiplicity_idx_2{i}=SC_idx_2';
    % SC_Persistent_Number_Non_Zero_Eigenvalue_idx_2{i}=SC_idx_2';
    % SC_Persistent_Quasi_Wiener_Index_idx_2{i}=SC_idx_2';
    % SC_Persistent_Spanning_Tree_Number_idx_2{i}=SC_idx_2';
    % SC_Persistent_Standard_Deviation_idx_2{i}=SC_idx_2';
    SC_Persistent_Generalized_Mean_Graph_Energy_Centroids_2{i}=SC_Centroids_2;
    % SC_Persistent_Laplacian_Graph_Energy_Centroids_2{i}=SC_Centroids_2;
    % SC_Persistent_Maximum_Centroids_2{i}=SC_Centroids_2;
    % SC_Persistent_Mean_Centroids_2{i}=SC_Centroids_2;
    % SC_Persistent_Minimum_Centroids_2{i}=SC_Centroids_2;
    % SC_Persistent_Moment_Second_Order_Centroids_2{i}=SC_Centroids_2;
    % SC_Persistent_Multiplicity_Centroids_2{i}=SC_Centroids_2;
    % SC_Persistent_Number_Non_Zero_Eigenvalue_Centroids_2{i}=SC_Centroids_2;
    % SC_Persistent_Quasi_Wiener_Index_Centroids_2{i}=SC_Centroids_2;
    % SC_Persistent_Spanning_Tree_Number_Centroids_2{i}=SC_Centroids_2;
    % SC_Persistent_Standard_Deviation_Centroids_2{i}=SC_Centroids_2;
    subplot(2,1,2);
    scatter(SC_Y_t_SNE_2(SC_idx_2==1,1),SC_Y_t_SNE_2(SC_idx_2==1,2),'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','r','LineWidth',1);
    hold on
    scatter(SC_Y_t_SNE_2(SC_idx_2==2,1),SC_Y_t_SNE_2(SC_idx_2==2,2),'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','b','LineWidth',1);
    hold on
    scatter(SC_Y_t_SNE_2(SC_idx_2==3,1),SC_Y_t_SNE_2(SC_idx_2==3,2),'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','k','LineWidth',1);
    scatter(SC_Centroids_2(:,1),SC_Centroids_2(:,2),'kx');
    legend('Cluster 1','Cluster 2','Cluster 3','Centroids','Location','best','FontSize',12,'FontName','Times New Roman');
    hold off
    title('\it k\rm-means','FontSize',20,'FontName','Times New Roman');
    str1='D:\NTU\work\Literatures\Hi-C\Figure\t-SNE-assisted k-means\t-SNE_assisted_k-means_Persistent_Generalized_Mean_Graph_Energy_Chromosome_';
    % str1='D:\NTU\work\Literatures\Hi-C\Figure\t-SNE-assisted k-means\t-SNE_assisted_k-means_Persistent_Laplacian_Graph_Energy_Chromosome_';
    % str1='D:\NTU\work\Literatures\Hi-C\Figure\t-SNE-assisted k-means\t-SNE_assisted_k-means_Persistent_Maximum_Chromosome_';
    % str1='D:\NTU\work\Literatures\Hi-C\Figure\t-SNE-assisted k-means\t-SNE_assisted_k-means_Persistent_Mean_Chromosome_';
    % str1='D:\NTU\work\Literatures\Hi-C\Figure\t-SNE-assisted k-means\t-SNE_assisted_k-means_Persistent_Minimum_Chromosome_';
    % str1='D:\NTU\work\Literatures\Hi-C\Figure\t-SNE-assisted k-means\t-SNE_assisted_k-means_Persistent_Moment_Second_Order_Chromosome_';
    % str1='D:\NTU\work\Literatures\Hi-C\Figure\t-SNE-assisted k-means\t-SNE_assisted_k-means_Persistent_Multiplicity_Chromosome_';
    % str1='D:\NTU\work\Literatures\Hi-C\Figure\t-SNE-assisted k-means\t-SNE_assisted_k-means_Persistent_Number_Non_Zero_Eigenvalue_Chromosome_';
    % str1='D:\NTU\work\Literatures\Hi-C\Figure\t-SNE-assisted k-means\t-SNE_assisted_k-means_Persistent_Quasi_Wiener_Index_Chromosome_';
    % str1='D:\NTU\work\Literatures\Hi-C\Figure\t-SNE-assisted k-means\t-SNE_assisted_k-means_Persistent_Spanning_Tree_Number_Chromosome_';
    % str1='D:\NTU\work\Literatures\Hi-C\Figure\t-SNE-assisted k-means\t-SNE_assisted_k-means_Persistent_Standard_Deviation_Chromosome_';
    if(i<=22)
        str2=num2str(i);
    elseif(i==23)
        str2='X';
    elseif(i==24)
        str2='Y';
    end
    Path=[str1,str2];
    print(Path,'-dtiffn','-r300');
end