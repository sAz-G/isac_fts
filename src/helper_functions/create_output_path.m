function output_path = create_output_path(name_output_folder)
%CREATE_OUTPUT_PATH Summary of this function goes here
%   Detailed explanation goes here

    if ispc
        output_path = fullfile('..\..\..\figures', name_output_folder);
    elseif isunix
        output_path = fullfile('../../../figures', name_output_folder);
    else
        output_path = fullfile('..\..\..\figures', name_output_folder);
    end 
    
    if ~exist(output_path, 'dir')
        mkdir(output_path)
    end

end