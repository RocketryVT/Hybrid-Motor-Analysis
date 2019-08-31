function include_repo(folder_path)
% Adds to the MATLAB path all folders in the repo
% 
% INPUTS
% folder_path - string (optional)
%               Directory path to the repo folder
% 
% @author: Matt Marti
% @date: 2019-04-18

% If no argument, add folders from here
if nargin == 0
    folder_path = '.';
end

% Add '/' to the end of the path name
if folder_path(end) ~= filesep
    folder_path = [ folder_path, filesep ];
end

% Recursively add directories to the path
recursive_add( [folder_path, '.'] );

end


function recursive_add(str)
% Recursively add all directories in str to the path
addpath(str);
dirstruct = dir(str);
for i = 3:length(dirstruct)
    d = dirstruct(i);
    dirname = [str filesep d.name];
    if exist(dirname, 'dir')
        recursive_add( dirname );
    end
end

end




