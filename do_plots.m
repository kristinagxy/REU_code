function [dict] = do_plots(mypath)

% plot PrshowX
% mypath: directory name that stores the .mat files
% files in mypath example: tomo_single_64_default_0.01_info.mat

d = uigetdir(mypath);
fileList = dir(fullfile(d, '*.mat'));

count = 1;

for i=1:length(fileList)
    name = string(fileList(i).name);
    if contains(name,'x')
        xpath = string(fileList(i).folder)+'/'+name;
        a = load(xpath);
        infoname = fieldnames(a);
        X = a.(infoname{1});
        prob_path = string(fileList(i).folder)+'/'+strrep(name,'x','prob');
        a = load(prob_path);
        infoname = fieldnames(a);
        ProbInfo = a.(infoname{1});
        setting = split(name,'_');
        label = string(strjoin(setting(1:5)));
        figure(count), clf, axes('FontSize', 18), hold on;
        idx = size(X,2);
        while sum(isnan(X(:,idx))) ~= 0 %|| sum(isinf(X(:,idx))) ~= 0
            idx=idx-1;
        end
        PRshowx(X(:,idx), ProbInfo), hold on 
        title(label);
        count=count+1;
    end
end