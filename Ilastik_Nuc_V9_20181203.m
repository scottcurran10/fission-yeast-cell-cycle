%Name: Ilastik_Nuc_V9
%Date last modified: 2018/12/03, debugging
%Author: Gautam Dey
%Collaborators: Scott Curran
%Goal: Extract information from Ilastik segmentations; with nuclear mask; includes binning 

%%%%%%%Notes%%%%%%
%Note for future: combination of Canny and Prewitt edge detectors
%actually do a great job of picking out both internal and external
%features- in case we need to replace Ilastik with deterministic
%segmenting at any point

%initialise
clear; close all

%load files
namefrag='sc317_cdc25-mNG';
indir='/Users/Gautam/Dropbox/Collaborations/Pombe Segmentation files/Nuclei Segmentation/';

%filenames
fname_bf=[namefrag '_BF.tif'];
fname_ng=[namefrag '_MAX3-7.tif'];
fname_mask1=[namefrag '_seg.tif'];
fname_mask2=[namefrag '_nuc.tif']; %nuclear mask
fname_mask_output=[namefrag '_seg_split.tif'];

%number of frames in stack
numframes=27;

%analysis parameters
solidity_thresh=0.959; %for kinked cells
large_disk=20; small_disk=12; %septation filter
minsize=4000; %minimum cell size
pixelsize=0.0651;
bin_min=5; bin_max=25; bin_step=0.25;
background_value=0; %change if not wildtype
splitornot=0; %1 for split, 0 for not
botperc=2; topperc=98; %plotting percentiles 

%open file for writing 
fileID1=fopen([indir namefrag '_raw_data.txt'],'w'); %15 columns, currently 
fprintf(fileID1,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','Frame','Cell','Frame_Cell','Length','Area','Intensity','Mean','Num_Nuclei','Nuc_Area','Nuc_Int','Nuc_Mean','Cyto_Area','Cyto_Int','Cyto_Mean','Nuc/Total');

%loop through stack
totdata_cell=[]; totdata_nuc=[]; cntr1=0; 
for f=1:(numframes)
    
    %input
    
    %mask
    cellmask1=tiffread([indir fname_mask1],f);
    cellmask1=cellmask1.data; cellmask1=~cellmask1;
    
    %nuclear mask
    nucmask=tiffread([indir fname_mask2],f);
    nucmask=nucmask.data; nucmask=~nucmask;
    nucmasklbl=bwlabel(nucmask);
      
    %fluorescence
    cellng=tiffread([indir fname_ng],f);
    cellng=cellng.data;
    
    %background subtraction
    cellng_sub=cellng-background_value; %negative values eliminated
    
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
    
    %cytoplasm only
    cellmask1cyto=cellmasklbl;
    cellmask1cyto(nucmask>0)=0;
    
    if f==1
        imwrite(cellmasklbl,[indir fname_mask_output],'Compression','none');
    else
        imwrite(cellmasklbl,[indir fname_mask_output],'Compression','none','WriteMode','append');
    end
  
    %cell length + intensity measurements
    celldata=regionprops(cellmask1split,'Extrema','PixelIdxList','Area');
    
    %cytoplasm measurements
    cytodata=regionprops(cellmask1cyto,'PixelIdxList','Area');
    
    %nuclear measurements
    nucdata=regionprops(nucmask,'PixelIdxList','Area');
    
    for c=1:length(celldata)
        
        %identifying matching nuclei
        nucids=unique(nucmasklbl(cellmasklbl==c));
        nucids=nucids(nucids>0);
        
        %Finding longest axis
        cext=celldata(c).Extrema;
        D=pdist(cext); clength=max(D);
        
        %intensity in cell (nucleus + cytoplasm)
        cint=cellng_sub(celldata(c).PixelIdxList);
        cint=sort(cint,'Descend');
        
        %intensity in cytoplasm
        cytoint=sum(cellng_sub(cytodata(c).PixelIdxList));
        cytoarea=cytodata(c).Area;
        
        %looping through nuclei
        if ~isempty(nucids)
            
            %track
            cntr1=cntr1+1;
            
            nucint=[]; nucarea=[]; 
            for n=1:length(nucids)
                
                nucarea=[nucarea; nucdata(nucids(n)).Area];
                nucidx=nucdata(nucids(n)).PixelIdxList;
                nucint=[nucint; sum(cellng_sub(nucidx))];
            end
            nucint=sum(nucint); nucarea=sum(nucarea); %combining nuclei
            
            %storing values;
            %totdata_cell=[totdata_cell; 1:f 2:cntr1 3:c 4:clength.*pixelsize 5:celldata(c).Area 6:sum(cint) 7:mean(cint) 8:length(nucids) 9:nucarea 10:nucint 11:nucint./nucarea 12:cytoarea 13:cytoint 14:cytoint./cytoarea 15:nucint./sum(cint)];
            totdata_cell=[totdata_cell; f cntr1 c clength.*pixelsize celldata(c).Area sum(cint) mean(cint) length(nucids) nucarea nucint nucint./nucarea cytoarea cytoint cytoint./cytoarea nucint./sum(cint)]; %#ok<*AGROW>
            fprintf(fileID1,'%6.2f\t%6.2f\t%6.2f\t%6.2f\t%10.2f\t%12.3f\t%12.3f\t%4.0f\t%10.2f\t%12.3f\t%10.2f\t%10.2f\t%12.3f\t%10.2f%12.3f\n',totdata_cell(cntr1,:));
            disp([f c])
            
        end
    end
end
fclose(fileID1);


%%

%%%% plotting routines %%%%%

%total data binning
fileID2=fopen([indir namefrag '_binned_data_singlenuc.txt'],'w');
fprintf(fileID2,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','Bin_Centre','Bin_Count','Intensity','Avg_Intensity','Nuc_Intensity','Nuc_Mean_Intensity','Nuc/Total_Ratio','Intensity_SD','Avg_Intensity_SD','Nuc_Intensity_SD','Nuc_Mean_Intensity_SD','Nuc/Total_Ratio_SD');

fileID3=fopen([indir namefrag '_binned_data_multinuc.txt'],'w');
fprintf(fileID3,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','Bin_Centre','Bin_Count','Intensity','Avg_Intensity','Nuc_Intensity','Nuc_Mean_Intensity','Nuc/Total_Ratio','Intensity_SD','Avg_Intensity_SD','Nuc_Intensity_SD','Nuc_Mean_Intensity_SD','Nuc/Total_Ratio_SD');

%plot vectors
%order: [1:f 2:cntr1 3:c 4:clength.*pixelsize 5:celldata(c).Area 6:sum(cint)
%7:mean(cint) 8:length(nucids) 9:nucarea 10:nucint 11:nucint./nucarea 12:cytoarea 13:cytoint 14:cytoint./cytoarea 15:nucint./sum(cint)]
singlenuc=find(totdata_cell(:,8)==1);
multinuc=find(totdata_cell(:,8)>1);

lvec1=totdata_cell(singlenuc,4);
ivec1=totdata_cell(singlenuc,6);
ivec_mn1=totdata_cell(singlenuc,7);
ivec_nuc1=totdata_cell(singlenuc,10); 
ivec_nuc_mn1=totdata_cell(singlenuc,11);
ivec_ratio1=totdata_cell(singlenuc,15);

lvec2=totdata_cell(multinuc,4);
ivec2=totdata_cell(multinuc,6);
ivec_mn2=totdata_cell(multinuc,7);
ivec_nuc2=totdata_cell(multinuc,10); 
ivec_nuc_mn2=totdata_cell(multinuc,11);
ivec_ratio2=totdata_cell(multinuc,15);

%binning by length
%binvec=(1.1*min(lvec)):0.25:(0.9*max(lvec));
binvec=bin_min:bin_step:bin_max;
binmat_mn1=zeros(length(binvec)-1,5); binmat_sd1=zeros(length(binvec)-1,5); bcent=zeros(length(binvec)-1,1);
binmat_mn2=binmat_mn1; binmat_sd2=binmat_sd1; 

for b=1:(length(binvec)-1)
    
    %bin centres
    bcent(b)=binvec(b)+((binvec(b+1)-binvec(b))./2);
    
    %binning each vector
    bin_ind1=find(lvec1<binvec(b+1) & lvec1>=binvec(b));
    binmat_mn1(b,:)=[mean(ivec1(bin_ind1)) mean(ivec_mn1(bin_ind1)) mean(ivec_nuc1(bin_ind1)) mean(ivec_nuc_mn1(bin_ind1)) mean(ivec_ratio1(bin_ind1))];
    binmat_sd1(b,:)=[std(ivec1(bin_ind1)) std(ivec_mn1(bin_ind1)) std(ivec_nuc1(bin_ind1)) std(ivec_nuc_mn1(bin_ind1)) std(ivec_ratio1(bin_ind1))];
    
    bin_ind2=find(lvec2<binvec(b+1) & lvec2>=binvec(b));
    binmat_mn2(b,:)=[mean(ivec2(bin_ind2)) mean(ivec_mn2(bin_ind2)) mean(ivec_nuc2(bin_ind2)) mean(ivec_nuc_mn2(bin_ind2)) mean(ivec_ratio2(bin_ind2))];
    binmat_sd2(b,:)=[std(ivec2(bin_ind2)) std(ivec_mn2(bin_ind2)) std(ivec_nuc2(bin_ind2)) std(ivec_nuc_mn2(bin_ind2)) std(ivec_ratio2(bin_ind2))];
    
    
    %print to file
    fprintf(fileID2,'%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\n',[bcent(b) length(bin_ind1) binmat_mn1(b,:) binmat_sd1(b,:)]);
    fprintf(fileID3,'%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\n',[bcent(b) length(bin_ind2) binmat_mn2(b,:) binmat_sd2(b,:)]);
   
end
fclose(fileID2); 

%%

%colormaps
c1=([211 224 237])./255;
c2=([45 98 173])./255;
c3=([44 162 95])./255;
c4=([153 216 201])./255;

%set x axis
xmin=bin_min; xmax=bin_max;

%plotting
figure;
subplot(5,1,1); hold on
plot(lvec1,ivec1,'.','MarkerFaceColor',c1,'MarkerEdgeColor',c1); errorbar(bcent,binmat_mn1(:,1),binmat_sd1(:,1),'.','Color',c2,'MarkerFaceColor',c2,'MarkerEdgeColor',c2);
plot(lvec2,ivec2,'.','MarkerFaceColor',c3,'MarkerEdgeColor',c3); errorbar(bcent,binmat_mn2(:,1),binmat_sd2(:,1),'.','Color',c4,'MarkerFaceColor',c4,'MarkerEdgeColor',c4);
axis([xmin xmax prctile(ivec1,botperc) prctile(ivec1,topperc)]);
title(namefrag); xlabel('Cell length'); ylabel('Intensity');

subplot(5,1,2); hold on
plot(lvec1,ivec_mn1,'.','MarkerFaceColor',c1,'MarkerEdgeColor',c1); errorbar(bcent,binmat_mn1(:,2),binmat_sd1(:,2),'.','Color',c2,'MarkerFaceColor',c2,'MarkerEdgeColor',c2);
plot(lvec2,ivec_mn2,'.','MarkerFaceColor',c3,'MarkerEdgeColor',c3); errorbar(bcent,binmat_mn2(:,2),binmat_sd2(:,2),'.','Color',c4,'MarkerFaceColor',c4,'MarkerEdgeColor',c4);
axis([xmin xmax prctile(ivec_mn1,botperc) prctile(ivec_mn1,topperc)]);
title(namefrag); xlabel('Cell length'); ylabel('Average intensity');

subplot(5,1,3); hold on
plot(lvec1,ivec_nuc1,'.','MarkerFaceColor',c1,'MarkerEdgeColor',c1); errorbar(bcent,binmat_mn1(:,3),binmat_sd1(:,3),'.','Color',c2,'MarkerFaceColor',c2,'MarkerEdgeColor',c2);
plot(lvec2,ivec_nuc2,'.','MarkerFaceColor',c3,'MarkerEdgeColor',c3); errorbar(bcent,binmat_mn2(:,3),binmat_sd2(:,3),'.','Color',c4,'MarkerFaceColor',c4,'MarkerEdgeColor',c4);
axis([xmin xmax prctile(ivec_nuc1,botperc) prctile(ivec_nuc1,topperc)]);
title(namefrag); xlabel('Cell length'); ylabel('Nuclear Intensity');

subplot(5,1,4); hold on
plot(lvec1,ivec_nuc_mn1,'.','MarkerFaceColor',c1,'MarkerEdgeColor',c1); errorbar(bcent,binmat_mn1(:,4),binmat_sd1(:,4),'.','Color',c2,'MarkerFaceColor',c2,'MarkerEdgeColor',c2);
plot(lvec2,ivec_nuc_mn2,'.','MarkerFaceColor',c3,'MarkerEdgeColor',c3); errorbar(bcent,binmat_mn2(:,4),binmat_sd2(:,4),'.','Color',c4,'MarkerFaceColor',c4,'MarkerEdgeColor',c4);
axis([xmin xmax prctile(ivec_nuc_mn1,botperc) prctile(ivec_nuc_mn1,topperc)]);
title(namefrag); xlabel('Cell length'); ylabel('Avg Nuclear Intensity');

subplot(5,1,5); hold on
plot(lvec1,ivec_ratio1,'.','MarkerFaceColor',c1,'MarkerEdgeColor',c1); errorbar(bcent,binmat_mn1(:,5),binmat_sd1(:,5),'.','Color',c2,'MarkerFaceColor',c2,'MarkerEdgeColor',c2);
plot(lvec2,ivec_ratio2,'.','MarkerFaceColor',c3,'MarkerEdgeColor',c3); errorbar(bcent,binmat_mn2(:,5),binmat_sd2(:,5),'.','Color',c4,'MarkerFaceColor',c4,'MarkerEdgeColor',c4);
axis([xmin xmax prctile(ivec_ratio1,botperc) prctile(ivec_ratio1,topperc)]);
title(namefrag); xlabel('Cell length'); ylabel('Nuclear/Total Ratio');

%print figure to folder
print([indir namefrag '_plots'],'-dpdf','-r300','-fillpage')


