function [dict] = make_plots(mypath,type)

% plot Enorm
% mypath: directory name that stores the .mat files
% files in mypath example: tomo_single_64_default_0.01_info.mat
% type: 'precision','size','blurlevel','noise'

k=5;
if strcmp(type,'precision')
    k = 2;
elseif strcmp(type,'size')
    k = 3;
else strcmp(type,'blurlevel')
    k = 4;
end

d = uigetdir(mypath);
fileList = dir(fullfile(d, '*.mat'));
dict = containers.Map;
for i=1:length(fileList)
    name = string(fileList(i).name);
    if contains(name,'info')
        path = string(fileList(i).folder)+'/'+name;
        a = load(path);
        setting = split(name,'_');
        key = string(strjoin(setting(setdiff(1:5,k))));
        label = string(setting(k));
        if isKey(dict,key)
            record = dict(key);
            infoname = fieldnames(a);
            a = a.(infoname{1});
            record(label)=a.Enrm;
            dict(key)=record;
        else
            record = containers.Map;
            infoname = fieldnames(a);
            a = a.(infoname{1});
            record(label)=a.Enrm;
            dict(key)=record;
        end
    end
end
n=1;
thekeys = keys(dict);
shape = {'-o','-d','-<'};
for i = 1:length(thekeys)
    thetitle = thekeys{i};
    a = dict(thetitle);
    figure(n), clf, axes('FontSize', 18), hold on
    title(thetitle)
    leskeys = keys(a);
    for j=1:length(leskeys)
        plot(a(leskeys{j}),shape{j}, 'LineWidth',1,'DisplayName',leskeys{j}),hold on
        ylabel('Enrm')
        xlabel('iteration')
        legend();
    end
    n=n+1;
    
end