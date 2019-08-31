function Test_All()
% Test_All.m
% 
% Recursively searches the Matlab folder in the Hybrid-Motor-Analysis repo
% and runs scripts that begin with the characters "test_" (Case
% insensitive).
%
% @authors: Matt Marti
% @date: 2019-08-31

clear, clc, close all

%% Recursively search for test scripts

% Determine repo location in the directory tree
m = inmem('-completenames');
full_testall_path = m{contains(m, 'Test_All.m')};
idum = regexp(full_testall_path, filesep);
path = full_testall_path(1:idum(end));

% Recursively run test scripts
include_repo(path);
recursive_test(path);


end

function recursive_test(str)
% Recursively add all directories in str to the path

% Find and run tests


% Perform action in folders
dirstruct = dir(str);
for i = 3:length(dirstruct)
    d = dirstruct(i);
    dirname = [str filesep d.name];
    if exist(dirname, 'dir')
        recursive_test( dirname );
    end
end

end
