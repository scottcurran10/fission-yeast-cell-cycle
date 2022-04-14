%Name: ImageStream_Heatmaps_V4
%Date last modified: 2020/01/15
%Authors: Gautam Dey, Scott Curran
%Goal: Generate heatmaps from ImageStream data 

%%%%%%%

%initialise
clear; %close all

%load parameters
indir='/Users/gautam/Dropbox/Collaborations/Pombe Segmentation files/Heatmaps/Jan 2020/';

namefrag='S';
loadfulldata=0;
binfulldata=0;
sortdata=1;
numstrains=49;
numcutoff=500;
load ([indir 'StrainGroups.mat']);

%% concatenate
if loadfulldata==0
    
    load([indir 'ImageStream-RawData-20200115.mat']);
    %find extents
    FullStrainData=struct;
    for nStrain=1:numstrains
        
        StrainName=['= ' namefrag num2str(nStrain) ';'];
        TempVar=genvarname('TempVar');
        eval([TempVar StrainName]);
        
        FullStrainData(nStrain).data=TempVar;
        
        disp(nStrain)
    end
    
    save([indir 'FullStrainData.mat'],'FullStrainData');
end

%% Various binning

if binfulldata==0
    
    load([indir 'FullStrainData.mat']);
    
    StrainNameVec={}; LengthBinTot=[]; SeptVec=[];
    for nStrain=1:numel(FullStrainData)
        
        D=FullStrainData(nStrain).data;
        Dmat=D{:,1:4};
        SeptVec=[SeptVec; D{3,7}];
        SN=D.Properties.VariableNames;
        StrainNameVec{nStrain}=SN{7}; %#ok<*SAGROW>
        %keyboard
        
        LengthBinTot=[LengthBinTot; Dmat(:,1)]; %#ok<*AGROW>
        
    end
    LengthBinTot=unique(LengthBinTot);
    LengthBinTot=LengthBinTot(~isnan(LengthBinTot));
    
    %correction
    StrainNameVec{22}='SC265';
    
    %Concatenate all strains
    IntMat=[]; IntMatRaw=[]; SDMat=IntMat; CountMat=IntMat; StrCntr=0;
    for nStrain=1:numel(FullStrainData)
        
        %strain matching
        StrMatch=strcmp(Strain,StrainNameVec{nStrain});
        %keyboard
        if sum(StrMatch)==1
            
            StrMatch=find(StrMatch==1);
            StrCntr=StrCntr+1;
            
            D=FullStrainData(nStrain).data; Dmat=D{:,1:4};
            IntMatRaw(nStrain)=sum(Dmat(:,2));
            %keyboard
            for nBin=1:numel(LengthBinTot)
                
                binidx=find(Dmat(:,1)==LengthBinTot(nBin));
                
                %get minimum value for normalization
                minidx=find(Dmat(:,4)>500);
                minval=min(Dmat(minidx,2));
                
                if ~isempty(binidx) && Dmat(binidx,4)>500 %remove bins <500
                    
                    IntMat(StrMatch,nBin)=Dmat(binidx,2)./minval;
                    SDMat(StrMatch,nBin)=Dmat(binidx,3)./minval;
                    CountMat(StrMatch,nBin)=Dmat(binidx,4)./sum(Dmat(:,4));
                    
                else
                    
                    IntMat(StrMatch,nBin)=NaN;
                    SDMat(StrMatch,nBin)=NaN;
                    CountMat(StrMatch,nBin)=NaN;
                    
                end
            end
        end
    end
    save([indir 'RawHeatMap.mat'],'IntMat','SDMat','CountMat','LengthBinTot','StrainNameVec','SeptVec','IntMatRaw');
end
% subplot(1,2,1)
% imagesc(IntMat); colormap gray
% subplot(1,2,2)
% imagesc(CountMat); colormap gray


%%

if sortdata==0
    
    load([indir 'RawHeatMap.mat']);
    
    %sort by count
    %IntCountMat=zeros(numstrains,60);
    IntSeptMat=zeros(numstrains,60);
    for nStrain=1:numstrains
        
        [mval,countidx]=max(CountMat(nStrain,:));
        septidx=find(LengthBinTot>SeptVec(nStrain),1);
        
        %         Stval1=30-countidx; Enval1=(Stval1+numel(LengthBinTot))-1;
        %         IntCountMat(nStrain,Stval1:Enval1)=IntMat(nStrain,:);
        %
        Stval2=50-septidx; Enval2=(Stval2+numel(LengthBinTot))-1;
        IntSeptMat(nStrain,Stval2:Enval2)=IntMat(nStrain,:);
        
        %sortvec=1:size(CountMat,2); sortvec=sortvec-countidx;
        %plot(sortvec,CountMat(nStrain,:),'-'); hold on
        %plot(sortvec,IntMat(nStrain,:),'-'); hold on
        
        
    end
    
end

%% plotting routines

%sets (see notes)
keepset1=find(Group==1);
keepset2=find(Group==2); keepset2=setdiff(keepset2,[34 38]);
keepset3=find(Group==3);

%keepset1=setdiff(1:47,[34 38]);

figure;
imagesc(IntMat(keepset1,2:49)); colorbar; title ('Aligned by length'); caxis([0.8 2]); %caxis([0.8 2.5]);
yticks(1:length(keepset1));
yticklabels(Genotype(keepset1));
xticktotalvec=LengthBinTot(2:49); %xticktotalvec=xticktotalvec+1; 
xticks([2 5 8 11 14 17 20 23 26 29 32 35 38 41 44 47]+1);
xticklabels(xticktotalvec(xticks));

%colormap(brewermap([],'*RdBu'))
colormap(brewermap([],'*Spectral'))
%colormap(brewermap([],'*PuBu'))

figure; subplot(2,1,1)
imagesc(IntMat(keepset2,2:49)); colorbar; title ('Aligned by length'); caxis([0.8 2]); %caxis([0.8 2.5]);
yticks(1:length(keepset2));
yticklabels(Genotype(keepset2));
xticktotalvec=LengthBinTot(2:49); %xticktotalvec=xticktotalvec+1; 
xticks([2 5 8 11 14 17 20 23 26 29 32 35 38 41 44 47]+1);
xticklabels(xticktotalvec(xticks));

%colormap(brewermap([],'*RdBu'))
%colormap(brewermap([],'*PuBu'))
colormap(brewermap([],'*Spectral'))

subplot(2,1,2)
imagesc(IntMat(keepset3,2:49)); colorbar; title ('Aligned by length'); caxis([0.8 2.5]); %caxis([0.8 2.5]);
yticks(1:length(keepset3));
yticklabels(Genotype(keepset3));
xticktotalvec=LengthBinTot(2:49); %xticktotalvec=xticktotalvec+1; 
xticks([2 5 8 11 14 17 20 23 26 29 32 35 38 41 44 47]+1); 
xticklabels(xticktotalvec(xticks));

%colormap(brewermap([],'*RdBu'))
%colormap(brewermap([],'*PuBu'))
colormap(brewermap([],'*Spectral'))
