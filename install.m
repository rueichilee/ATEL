%% ATEL Package Installer
% Description: Adds the ATEL package and its submodules to your MATLAB path.

function install()
    fprintf('Installing ATEL Toolbox...\n');

    % 1. Locate the Package Root
    current_file = mfilename('fullpath');
    [root_path, ~, ~] = fileparts(current_file);
    
    % 2. Define Paths
    src_path = fullfile(root_path, 'src');
    exa_path = fullfile(root_path, 'example'); 
    test_path = fullfile(root_path, 'test');

    % 3. Check and Add 'src' 
    if exist(src_path, 'dir')
        addpath(genpath(src_path));
        fprintf('[OK] Added source modules: %s\n', src_path);
    else
        warning('Folder "src" not found at %s', src_path);
    end

    % 4. Check and Add 'example' 
    % FIX: Use the variable 'exa_path' we defined above
    if exist(exa_path, 'dir')
        addpath(genpath(exa_path));
        fprintf('[OK] Added example/data:  %s\n', exa_path);
    else
        % Fallback: Try looking for "examples" (plural) just in case
        exa_path_alt = fullfile(root_path, 'examples');
        if exist(exa_path_alt, 'dir')
            addpath(genpath(exa_path_alt));
            fprintf('[OK] Added examples/data:  %s\n', exa_path_alt);
        else
            warning('Example folder not found. Data may not load automatically.');
        end
    end
    
    % 5. Check and Add 'tests' 
    if exist(test_path, 'dir')
        addpath(genpath(test_path));
        fprintf('[OK] Added test:          %s\n', test_path);
    else
        warning('Tests folder not found at %s', test_path);
    end

    % 6. Add Root Path
    addpath(root_path);


    fprintf('-----------------------------------------\n');
    fprintf('ATEL Toolbox installed successfully!\n');
    fprintf('Run "test_simulation" to verify installation.\n');
    fprintf('-----------------------------------------\n');
end