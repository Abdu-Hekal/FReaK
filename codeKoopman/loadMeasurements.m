function [X,U,time] = loadMeasurements(path)
% load measurements and bring them to the correct format

    content = dir(path);
    X = {}; U = {}; time = {};

    for i = 1:length(content)
       if content(i).isdir && startsWith(content(i).name,'measurement')
           newPath = fullfile(path,content(i).name);
           input = csvread(fullfile(newPath,'input.csv'));
           traj = csvread(fullfile(newPath,'trajectory.csv'));
           t = csvread(fullfile(newPath,'time.csv'));
           X{end+1} = traj'; U{end+1} = input'; time{end+1} = t;
       end
    end
end

