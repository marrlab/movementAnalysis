%run can take several minutes
%out folder is created in "/path/to/file/"
%input file1
file1='/path/to/file/FibroblastsD4.csv';
%plot displacement
[mov100_Fib_D4, trackMat]=calculate_movement(file1);
%plot the mean displacement over time on a grid
plot_Neighborhood_Similarity(trackMat,file1);
%input file2
file2='/path/to/file/AdipocytesD4.csv';
%plot displacement
[mov100_Ad_D4, trackMat]=calculate_movement(file2);
%plot splines comapring mean displacements of file1 and file2
plot_splines(mov100_Fib_D4,mov100_Ad_D4,file1)

%% calculate displacement
function [mov100, trackMat]=calculate_movement(file)

%read and prepare data
input=readtable(file);
[filepath,name,~] = fileparts(file);

% create output folder
if ~exist([filepath,'/out'], 'dir')
    mkdir([filepath,'/out'])
end

%shift trackIDs to start from 1
input.TrackID=input.TrackID-1000000000+1;
%get all unique trackIDs
ids = unique(input.TrackID);
%movement array to save displacements in
movements=cell(max([input.Time])-1,1);
%matrix to save tracks for later analysis
trackMat=[];
%iterate over all tracks
for t = 1: size(ids,1)
    %extract all time points of track t from the data
    track = find(input.TrackID==ids(t));
    % iterate of every time point j of track t
    for j = 1:numel(track)-1
        trackMat((t),input.Time(track(j)),:)=[input.PositionX(track(j)),input.PositionY(track(j)),input.PositionZ(track(j))];
        %check if tracks are found in censecutive time points
        if input.Time(track(j+1)) - input.Time(track(j)) ==1
            %calculate dispalcement between this timepoint and the next one
            displacement = sqrt((input.PositionX(track(j+1))-input.PositionX(track(j)))^2+...
                (input.PositionY(track(j+1))-input.PositionY(track(j)))^2+...
                (input.PositionZ(track(j+1))-input.PositionZ(track(j)))^2);
            movements{j}=[movements{j} displacement];
        end
    end
    
end

%% exclude wrong tracks
% keep only cells which show a movement below 100 microns. Higher values
% are not observed in these cells and can be excluded as track mismatching
mov100=cell(max([input.Time])-1,1);
for t = 1: numel(movements)
    m= [movements{t}];
    %divide by 15 to get mean movement per minute / velocity
    mov100{t} =m(m<100)/15;
end
%% plot
%plot Velocity over time
figure
distributionPlot((mov100(1:numel(mov100))),'histOpt',2,'showMM',2,'colormap',1-gray(64));
xticks(0:4:numel(mov100))
xticklabels(0:numel(mov100)/4)
ylim([0,2])
xlabel('Time (h)')
ylabel('Velocity (\mum/min)')
set(gca,'FontSize',16)
%save plot
saveas(gca,sprintf('%s_velocity.tif',[filepath,'/out/',name]));

%plot zoom in
figure
distributionPlot((mov100(16:36)),'histOpt',2,'showMM',2,'colormap',1-gray(64));
xticks(1:4:25)
xticklabels(4:9)
ylim([0,0.6])
xlabel('Time (h)')
ylabel('Velocity (\mum/min)')
set(gca,'FontSize',16)
%save plot
saveas(gca,sprintf('%s_velocity_zoom.tif',[filepath,'/out/',name]));

end

%% compare movement similarity
function plot_Neighborhood_Similarity(trackMat, file)
input=readtable(file);
%shift trackIDs to start from 1
input.TrackID=input.TrackID-1000000000+1;
%get all unique trackIDs
ids = unique(input.TrackID);
allSimilarities =[];
for t =1:size(trackMat,2)-1
    vectorSimilarity =[];
    tp1 = squeeze(trackMat(:,t,:));
    tp2 = squeeze(trackMat(:,t+1,:));
    %keep only cells appearing at both tps
    int12=intersect(input.TrackID(input.Time==t),input.TrackID(input.Time==t+1));
    [~,ind]=intersect(ids,int12);
    tp1=tp1(ind,:);
    tp2=tp2(ind,:);
    %create a neighborhood network
    tri = delaunay(tp1);
    Dnet = distMatOnNet(tri);
    %calculate angle between movement vectors for each cell j
    for j = 1:size(Dnet,1)
        %find all direct neighbors
        neighbors = find(Dnet(j,:) == 1);
        angles=[];
        % calculate the 3D angles between the movement of one cell and all
        % its neighbors respectively
        for k =1:numel(neighbors)
            a= tp2(j,:)-tp1(j,:);
            b= tp2(neighbors(k),:)-tp1(neighbors(k),:);
            angles = [angles,atan2d(norm(cross(a,b)), dot(a,b))];
        end
        %save the mean neighbor angle of cell j
        vectorSimilarity=[vectorSimilarity;[tp1(j,:),mean(angles)]];
    end
    %map the vector similarity back to the cells x and y coordinates
    [xi, yi] = meshgrid(linspace(min(vectorSimilarity(:,1)),max(vectorSimilarity(:,1))),linspace(min(vectorSimilarity(:,2)),max(vectorSimilarity(:,2))));
    zi = griddata(vectorSimilarity(:,1),vectorSimilarity(:,2),180-vectorSimilarity(:,4), xi,yi);
    allSimilarities(t,:,:)=zi;
    
end

%% plot
%plot the mean vector similarity to the coordinate system
figure('Position',[100, 100, 900, 900]);
colormap 'jet'
s=surf(xi,yi,squeeze(nanmean(allSimilarities)));
xlabel('X (\mum)')
ylabel('Y (\mum)')
caxis([90,180]);
view(0,90)
c = colorbar;
ylabel(c, 'Mean similarity (°)','FontSize',24)
set(gca,'FontSize',20)
%save
[filepath,name,~] = fileparts(file);
saveas(gca,sprintf('%s_neighborSimilarity.tif',[filepath,'/out/',name]));
end


%% plot splines
function plot_splines(mov100_1, mov100_2,file)
%prepare figure
figure;
hold on;

% average movement and plot
mov=cellfun(@mean,(mov100_1));
xx = 0:.25:numel(mov100_1);
yy = spline(1:numel(mov100_1),mov,xx);
h(1)=plot(xx,yy,'LineWidth',2,'DisplayName',"Fibroblasts day 4");

% average movement and plot
mov=cellfun(@mean,(mov100_2));
xx = 0:.25:numel(mov100_2);
yy = spline(1:numel(mov100_2),mov,xx);
h(2)=plot(xx,yy,'LineWidth',2,'DisplayName','Adipocytes day 4');

%modify axes
ylim([0,0.5])
xlim([0,48])
legend(h([2,1]))
xticks(0:4:numel(mov100_2))
xticklabels(0:numel(mov100_2)/4)
xlabel('Time (h)')
ylabel('Velocity (\mum/min)')
set(gca,'FontSize',16)
%save
[filepath,~,~] = fileparts(file);
saveas(gca,sprintf('%sSplines.tif',[filepath,'/out/']));
end

%% calculate distance matrix for a networks based on triangulation
function Dnet= distMatOnNet(tri)
net = zeros(max(max(tri)));
for i = 1:size(tri,1)
    net(tri(i,1),tri(i,2)) = 1;
    net(tri(i,2),tri(i,1)) = 1;
    net(tri(i,1),tri(i,3)) = 1;
    net(tri(i,3),tri(i,1)) = 1;
    net(tri(i,2),tri(i,3)) = 1;
    net(tri(i,3),tri(i,2)) = 1;
end
Dnet = all_shortest_paths(sparse( net));
end