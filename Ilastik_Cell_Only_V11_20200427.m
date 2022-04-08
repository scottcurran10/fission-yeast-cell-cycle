%Name: Ilastik_Cell_Only_V7
%Date last modified: 2020/04/07
%Author: Gautam Dey
%Collaborators: Scott Curran
%Goal: Extract information from Ilastik segmentations; no nuclear mask
%Latest change: Change top pixels to percentage of top pixels 


%%%%%%%Notes%%%%%%
%Note for future: combination of Canny and Prewitt edge detectors
%actually do a great job of picking out both internal and external
%features- in case we need to replace Ilastik with deterministic
%segmenting at any point

%init
clear; close all

%load files
namefrag='sc571 cdc22-M45 cdc25-mNG 36';
indir='/Volumes/SCOTT2/Cdc25 PAPER 2022/cdc25-mNG cdc22-M45 200626/sc571 cdc22-M45 cdc25-mNG/36/MATLAB/';

%filenames
fname_bf=[namefrag '_BF.tif'];
fname_ng=[namefrag '_MAX2-6.tif'];
fname_mask1=[namefrag '_NNseg.tif'];
fname_mask_output=[namefrag '_seg_split.tif'];

%stack size
numframes=28;

%analysis params
solidity_thresh=0.959;
large_disk=20;
small_disk=12;
minsize=1000;
pixelsize=0.0651;
bin_min=0; bin_max=40; bin_step=0.25;
background_value= 631.6233824; %change if not wildtype
splitornot=0; %1 for split, 0 for not
toppixelpercent=15; %approximating size of nucleus
botperc=0.5; topperc=99.5; %plotting percentiles 

%open file for writing
fileID1=fopen([indir namefrag '_raw_data.txt'],'w');
fprintf(fileID1,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','Frame','Cell','Frame_Cell','Length','Area','Intensity','Mean','Ttop_Pix','Ratio_Top_to_Total');

%loop through stack
totdata_cell=[]; totdata_nuc=[]; cntr1=0; cntr2=0;
for f=1:(numframes)
    
    %input
    
    %mask
    cellmask1=tiffread([indir fname_mask1],f);
    cellmask1=cellmask1.data; cellmask1=~cellmask1;
    
    %fluorescence
    cellng=tiffread([indir fname_ng],f);
    cellng=cellng.data;
    
    %background subtraction
    cellng_sub=cellng-background_value; %added conversion to double for negative values
    
    %Remove edge objects 2018/07/11
    cellmask1=imclearborder(cellmask1);
    
    %thresholding
    
    if splitornot == 1
        
        %First identify kinked cells
        kinkedcells=bwpropfilt(cellmask1,'Solidity',[0 solidity_thresh]);
        kinkmarker1=imerode(kinkedcells,strel('disk',large_disk));
        cellmarker1=cellmask1; cellmarker1(kinkedcells==1)=0;
        
        %marker-based watershed
        cellmarker1=imerode(cellmarker1,strel('disk',small_disk));
        cellmarker1=im2bw(cellmarker1+kinkmarker1);  %#ok<*IM2BW>
        cellwat1=watershed(bwdist(cellmarker1));
        cellmask1split=cellmask1; cellmask1split(cellwat1==0)=0;
        
        %clear edges and any remaining noise
        cellmask1split=bwareaopen(cellmask1split,minsize);
        cellmask1split=imclearborder(cellmask1split);
        
    else
        cellmask1split=cellmask1;
    end
    
    %save masks to multipage tiff
    cellmasklbl=bwlabel(cellmask1split); %2018/07/11
    cellmasklbl=uint16(cellmasklbl);
    
    
    if f==1
        imwrite(cellmasklbl,[indir fname_mask_output],'Compression','none');
    else
        imwrite(cellmasklbl,[indir fname_mask_output],'Compression','none','WriteMode','append');
    end
    
    %cell length + intensity measurements
    celldata=regionprops(cellmask1split,'Extrema','PixelIdxList','Area');
    
    for c=1:length(celldata)
        
        %track
        cntr1=cntr1+1;
        
        %Finding longest axis
        cext=celldata(c).Extrema;
        D=pdist(cext); clength=max(D);
        
        %intensity in cell (nucleus + cytoplasm)
        cint=cellng_sub(celldata(c).PixelIdxList);
        cint=sort(cint,'Descend');
        
        %Added for percentage top pixels instead of absolute number 
        tpp=toppixelpercent./100;
        numtoppixels=floor(tpp.*length(cint)); %keyboard
       
        %order: [frame_number total_count frame_count clength cell_area total_int mean_int top_pix ratio];
        totdata_cell=[totdata_cell; f cntr1 c clength.*pixelsize celldata(c).Area sum(cint) mean(cint) mean(cint(1:numtoppixels)) (mean(cint(1:numtoppixels)))./mean(cint)]; %#ok<*AGROW>
        disp([f c])
        fprintf(fileID1,'%6.2f\t%6.2f\t%6.2f\t%6.2f\t%10.2f\t%12.2f\t%12.3f\t%6.2f\t%6.2f\n',totdata_cell(cntr1,:));
        
        %keyboard
    end
end
fclose(fileID1);

% %Write to Excel- not working with MacOS & current MATLAB release
% xlswrite([indir namefrag 'data.xls'],{'Cell length','Area','Intensity','Top 20 pix','Top 5 pix'},1,'A1');
% xlswrite([indir namefrag 'data.xls'],totdata,1,'A2');

%%

%%%% plotting routines %%%%%

%total data binning
fileID2=fopen([indir namefrag '_binned_data.txt'],'w');
fprintf(fileID2,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','Bin_Centre','Bin_Count','Intensity','Avg_Intensity','Top','Top_ratio','Intensity_SD','Avg_Intensity_SD','Top_SD','Top_ratio_SD');
%order: [frame_number total_count frame_count clength cell_area total_int mean_int top20pix top5pix];

%plot vectors
lvec=totdata_cell(:,4);
ivec=totdata_cell(:,6);
ivec_mn=totdata_cell(:,7);
ivec_top=totdata_cell(:,8);
ivec_ratio=totdata_cell(:,9);

%binning by length
%binvec=(1.1*min(lvec)):0.25:(0.9*max(lvec));
binvec=bin_min:bin_step:bin_max;
binmat_mn=zeros(length(binvec)-1,4); binmat_sd=zeros(length(binvec)-1,4); bcent=zeros(length(binvec)-1,1);

for b=1:(length(binvec)-1)
    
    %bin centres
    bcent(b)=binvec(b)+((binvec(b+1)-binvec(b))./2);
    
    %binning each vector
    bin_ind=find(lvec<binvec(b+1) & lvec>=binvec(b));
    binmat_mn(b,:)=[mean(ivec(bin_ind)) mean(ivec_mn(bin_ind)) mean(ivec_top(bin_ind)) mean(ivec_ratio(bin_ind))];
    binmat_sd(b,:)=[std(ivec(bin_ind)) std(ivec_mn(bin_ind)) std(ivec_top(bin_ind)) std(ivec_ratio(bin_ind))];
    
    %print to file
    fprintf(fileID2,'%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\n',[bcent(b) length(bin_ind) binmat_mn(b,:) binmat_sd(b,:)]);
    
end
fclose(fileID2); 


%%

%colormaps
c1=([211 224 237])./255;
c2=([45 98 173])./255;

%set x axis
xmin=bin_min; xmax=bin_max;

%plotting
figure;
subplot(4,1,1); hold on
plot(lvec,ivec,'.','MarkerFaceColor',c1,'MarkerEdgeColor',c1); errorbar(bcent,binmat_mn(:,1),binmat_sd(:,1),'.','Color',c2,'MarkerFaceColor',c2,'MarkerEdgeColor',c2);
%axis([xmin xmax min(ivec) max(ivec)]);
axis([xmin xmax prctile(ivec,botperc) prctile(ivec,topperc)]);
title(namefrag); xlabel('Cell length'); ylabel('Intensity');

subplot(4,1,2); hold on
plot(lvec,ivec_mn,'.','MarkerFaceColor',c1,'MarkerEdgeColor',c1); errorbar(bcent,binmat_mn(:,2),binmat_sd(:,2),'.','Color',c2,'MarkerFaceColor',c2,'MarkerEdgeColor',c2);
%axis([xmin xmax min(ivec_mn) max(ivec_mn)]);
axis([xmin xmax prctile(ivec_mn,botperc) prctile(ivec_mn,topperc)]);
title(namefrag); xlabel('Cell length'); ylabel('Average intensity');

subplot(4,1,3); hold on
plot(lvec,ivec_top,'.','MarkerFaceColor',c1,'MarkerEdgeColor',c1); errorbar(bcent,binmat_mn(:,3),binmat_sd(:,3),'.','Color',c2,'MarkerFaceColor',c2,'MarkerEdgeColor',c2);
%axis([xmin xmax min(ivec_20) max(ivec_20)]);
axis([xmin xmax prctile(ivec_top,botperc) prctile(ivec_top,topperc)]);
title(namefrag); xlabel('Cell length'); ylabel('Top pixels');

subplot(4,1,4); hold on
plot(lvec,ivec_ratio,'.','MarkerFaceColor',c1,'MarkerEdgeColor',c1); errorbar(bcent,binmat_mn(:,4),binmat_sd(:,4),'.','Color',c2,'MarkerFaceColor',c2,'MarkerEdgeColor',c2);
%axis([xmin xmax min(ivec_5) max(ivec_5)]);
axis([xmin xmax prctile(ivec_ratio,botperc) prctile(ivec_ratio,topperc)]);
title(namefrag); xlabel('Cell length'); ylabel('Ratio top/total');

%print figure to folder
print([indir namefrag '_plots'],'-dpdf','-r300','-fillpage')


