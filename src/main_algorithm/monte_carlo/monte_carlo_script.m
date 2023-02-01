%%
clear;
clc;
close all;

addpath(genpath("..\..\..\src"));


%% Simulation parameter

run("call_hyperParam.m")

monte = 1; 
results_monte = [];

for k=1:3
     results_monte.(strcat('sim_loop',num2str(k))) = NaN;
end

%% save workspace
for mn = 1:monte
  results_monte.(strcat('sim_loop',num2str(mn))) = multi_stage(params,setup,results);  
end

%%
sv = 1;
folder_name = "default";
if sv == 1
    
    pth =  "..\..\..\results\";
    sfx = uint16(exist(pth+folder_name, 'dir')/7);
    
    if ~sfx
        mkdir(pth+folder_name);
        folder_name = pth+folder_name;
    else
        name = pth+folder_name;
        while exist(name, 'dir')
            name = pth+folder_name + string(sfx);
            if ~exist(name, 'dir')
                mkdir(name);
                folder_name = name;
                break;
            end
            sfx = sfx +1;
        end
    end
end


save(folder_name + '\workspace.mat')