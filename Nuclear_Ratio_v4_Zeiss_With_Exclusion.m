
%%=========================================================================
%
% Nuclear_Mask:
%          Segmentation of images based on DAPI signal to quantify
%          nuclear orienation, nuclear aspect ratio, nuclear area,
%          and average nuclear staining intensity for SIGNAL image.
%          Also calculates ratio of nuclear signal to cytoplasmic signal just
%          outside the nucleus (Ring of set diameter around nucleus).
%
% Inputs:         
%          DAPI images (16 bit tifs with format DAPI0001.tif)
%          SIGNAL images (16 bit tifs with format SIGNAL0001.tif)
%          ShadeCorrect.tif (16 bit tif of dye on shade channel, in base directory)  
%          imageno  (number of images per group)
%
% Editable Hardcoded Parameters:
%          thresh (threshold used for DAPI segmentation, 0-1, default 0.4)
%          minNucSize (minimum pixel size for each nucleus, default 50)
%          maxNucSize (maximum pixel size for each nucleus, default 1000)
%
% Outputs:
%          Data01.csv (csv file for each image)
%               Column 1:  Nuclear Aspect Ratio
%               Column 2:  Nuclear Orientation, 0-90 degrees, 90 = Vertical
%               Column 3:  Nuclear Area in pixels
%               Column 4:  X coordinate of nuclear centroid (0,0)= top left
%               Column 5:  Y coordinate of nuclear centroid
%               Column 6:  Average Signal in Nucleus (shade corrected)
%               Column 7:  Average Signal outside the nucleus
%               Column 8:  Ratio of nuclear to cytoplasmic signal
%
%          DAPI_Mask_0001.tif (numbered image of each masked nuclus)
%
%  Tristan Driscoll and Garrett McDaniel
%  tristan.driscoll@gmail.com
%  Driscoll Lab
%  8/24/23
%
%
%%
close all;
clear all;

imageno=20;     % # of images to analyze, format: DAPI0001 and YAP0001
ScaleFactor=0.5;  % To re-scale dim images, set to 1 for no rescale

%===================Initializations===================================
cmap=colormap('jet');  %set up color map
cmap(1,3)=0.0;  %set zero to black on colormap
warning('off','images:initSize:adjustingMag');  %turn off stupid warnings
thresh=[0.1]; %determines Canny thresh for Nuclear Segmentation
minNucSize=750;  %Min nucleus size used
maxNucSize=7000;  %Max nucleus size used
avgbk=0;
RingWidth=5;  %Size of the ring for cytoplasm signal
NARmax=2.5;

SigMin=5;
SigMax=4000;

%=================Input Parameters
%=========================================
%Directory Base
img_direc_base = uigetdir;  %user selects base directory
cd(img_direc_base);

subdirec='\Field_';  %Name of image subdirecs
Preface='';
DAPInames='DAPI';  %Name of dapi images
SIGnames='GREEN';    %Name of signal images

%Gradient Corection Images
% gradRAW=imread(strcat(img_direc_base,'/ShadeCorrect.tif'));
% gradRAW=double(gradRAW)-avgbk;
% gradNORM=gradRAW./max(max(gradRAW));



for j=1:imageno
    ExcludeRegion=zeros(1,1);
    time=1;
    clear NAR adjAngle NucArea adjAngle a NucMasks SigSum MeanSig SigNormArea data_out CytoArea CytoNormArea SigNormArea NCRatio
    cellsDirec=strcat(img_direc_base,'/',subdirec,num2str(j,'%02d'));
    cd(cellsDirec);
    DapiImgName = strcat(DAPInames, num2str(time,'%04d'),'.tif');
    DapiImg=imread(DapiImgName);
    name_DMask=strcat('DAPI_Mask_',num2str(time,'%04d'),'.tif');
    
    SigImgName = strcat(SIGnames, num2str(time,'%04d'),'.tif');
    SignalImg=imread(SigImgName);
    SignalImgOrig=SignalImg;
    SignalImgOrigTemp=SignalImgOrig;
    SignalImg=double(SignalImg);
    SignalImg=SignalImg-avgbk;
    %SignalImg=SignalImg./gradNORM;
    
    %maxDapi=max(max(DapiImg));
    
    DapiImg=imadjust(DapiImg,[0,ScaleFactor]);
    
    if max(DapiImg)>256
        DapiImg=im2uint8(DapiImg);
        
    end
    
    orig=DapiImg;
    
    %DapiImg=imadjust(DapiImg);
    BWs = edge(DapiImg,'canny', thresh); % edge detection; value of thresh determines sensitivity
    se90 = strel('line', 2, 90);
    se0 = strel('line', 2, 0);
    BWsdil = imdilate(BWs, [se90 se0]); % the borders selected by edge detection are turned into circles 3 px in diameter
    BWdfill = imfill(BWsdil, 'holes'); % fills areas that are not connected to the background
    BWnobord = imclearborder(BWdfill, 4); % removes nuclei intersecting the image edge
    
    seD = strel('disk',2);
    BWfinal = imerode(BWnobord,seD); %erodes edges of shapes
    figTest=figure;
    imshow(BWfinal);
    
    
    
    % initializations
    k=1;
    a=[0,0];
    
    
    %=====================================================================================
    %=====================================================================================
    %
    %             Nuclear SEGMENTATION
    %
    %=====================================================================================
    %=====================================================================================
    tic
    % cluster each nucleus and determine centroid
    [nucleus_cluster_matrix, nucleus_cluster_number] = bwlabeln(BWfinal);
    
    % determines centroid of each nucleus, uses PCA to find the long and short
    % axis of the shape, and finds the angle between the long axis a
    excluding=1;
    while excluding
        for i=1:nucleus_cluster_number
            [row,col] = find(nucleus_cluster_matrix == i);
            if length(row)> minNucSize && length(row)< maxNucSize % exclude any cluster that's too small or big
                if sum(ExcludeRegion(1,:)==i)==0
                    nucleus_cluster_centroid{k} = [round(mean(row)),round(mean(col))];
                    nucleus_set{k} = [row,col];
                    total_matrix = zeros(2);
                    for r=1:length(nucleus_set{k})
                        point_vector = nucleus_set{k}(r,:)-nucleus_cluster_centroid{k};
                        point_matrix = point_vector.'*point_vector;
                        total_matrix = total_matrix + point_matrix;
                    end
                    % find eigenvalues/vectors and sort
                    [eigenvectors,eigenvalues] = eig(total_matrix);
                    if max(eigenvalues(:,1))>max(eigenvalues(:,2))
                        long_axis = max(eigenvalues(:,1))^.5;
                        short_axis = max(eigenvalues(:,2))^.5;
                    else
                        long_axis = max(eigenvalues(:,2))^.5;
                        short_axis = max(eigenvalues(:,1))^.5;
                    end
                    % NAR and angle calulation


                    NucArea(k) = length(row);
                    NAR(k) = long_axis/short_axis;
                    angle(k) = acosd(dot(eigenvectors(:,1),[1;0]));
                    adjAngle(k)=abs(abs(angle(k)-90)-90); %angle from Horizontal
                    a(k,:)=nucleus_cluster_centroid{k};



                    k=k+1;
                    clear total_matrix point_matrix point_vector long_axis short_axis
                else
                    
                    for l=1:length(row)
                        BWfinal(row(l),col(l))=0; % remember, matrix vs img coord
                    end
                end


            else
                for l=1:length(row)
                    BWfinal(row(l),col(l))=0; % remember, matrix vs img coord
                end
            end

            clear row col
        end

        mask=(DapiImg==255);
        mask1=(DapiImg==254);
        DapiImg(mask)=252;
        DapiImg(mask1)=252;
        BWoutline = bwperim(BWfinal);
        Segout = DapiImg; Segout(BWoutline) = 255; %perim outline on original image
        
        test=imoverlay(SignalImgOrig,BWoutline,[1,0,1]);
        figTest=figure;
        imshow(test)
        
        sprintf(['Click all regions you would like to exclude for image %i' ...
            '\nIf you would not like to exclude any region or are finished selecting, press enter.\n'...
            'If you are having difficult seeing regions to exclude, press I then enter to invert the image\n'],j)      
        press=0;   
        
       [x,y,button]=ginput;
       if ~isempty(x)&&sum(button~=105)==length(button)   
                for s=1:length(x)

                    if nucleus_cluster_matrix(floor(y(s,1)),floor(x(s,1)))~=0
                        newExcludeRegion=nucleus_cluster_matrix(floor(y(s,1)),floor(x(s,1)));
                        ExcludeRegion=[ExcludeRegion,newExcludeRegion];
                        
                    else
                        sprintf('Selection %i was an invalid region',s)
                    end  
                end
                clear nucleus_set NucArea NAR angle adjAngle a nucleus_cluster_centroid x y
                        k=1;
       elseif isempty(button)
           excluding=0; 
        end 
        
        if button==105
            SignalImgOrigTemp=max(max(SignalImgOrigTemp))-SignalImgOrigTemp;
            test=imoverlay(SignalImgOrigTemp,BWoutline,[1,0,1]);
            imshow(test)
            clear x y 
        end
            
         
        
    end
    
    NucArea(NucArea==0)=[];
    
    
    
    %refresh to eliminate filtered objects
    [nucleus_cluster_matrix_refresh, nucleus_cluster_number_refresh] = bwlabeln(BWfinal);
    
    % OUTPUTS
    centroids=a(:,:);
    
    %save original image with text overlay
    
    %         SegoutSignal=SignalImg;
    %         SegoutSignal(BWoutline)=0;
    
    
    fig1=figure;
    SignalImgScaled=imadjust(SignalImgOrig, [0 ScaleFactor], [0 1]);
    
    SegoutSignal=imoverlay(SignalImgScaled, BWoutline, [1 0 1]);
    imshow(SegoutSignal)
    hold on;
    
    imagesc(SegoutSignal, [SigMin SigMax])
    colormap('gray');
    
    title('outlined original YAP image');
    %title(DapiFileName);
    
    
    for li=1:length(NucArea)
        text(centroids(li,2) , centroids(li,1), num2str(li), 'Color', 'm');
        hold on
    end
    axis equal;
    axis off;
    

    saveas(fig1, name_DMask, 'tif')
    data_out = [NAR', adjAngle', NucArea',centroids];
    
    
    toc
    
    %=====================================================================================
    %=====================================================================================
    %
    %            Signal SEGMENTATION
    %
    %
    %=====================================================================================
    %=====================================================================================
    
    % NUCLEAR SIGNAL
    tic
    fig2=figure;
    AvgSignal=mean(mean(SignalImg));
    MaxSignal=max(max(SignalImg));
    %imshow(SignalImg);
    colormap(cmap);
    colorbar;
    
    origSig=SignalImg;
    [xSig,ySig]=size(SignalImg);
    
    SigSum=zeros(nucleus_cluster_number_refresh, 1);
    BackSum=zeros(nucleus_cluster_number_refresh, 1);
    
    MeanSig=mean(SignalImg);
    
    for p=1:nucleus_cluster_number_refresh
%         if sum(selectImageRegion(j,2:end)==p)==0
            for n=1:xSig
                for m=1:ySig
                    if nucleus_cluster_matrix_refresh(n,m)==p
                        SigSum(p)=SigSum(p)+SignalImg(n,m);
                    end
                end
            end
%         end
    end

    SigNormArea=SigSum./NucArea';
    
    
    
    %Make Nuclear Ring Masks
    
    NucMasks=zeros(xSig,ySig,nucleus_cluster_number_refresh);
    RingMasks=zeros(xSig,ySig,nucleus_cluster_number_refresh);
    
    for p=1:nucleus_cluster_number_refresh
%         if sum(selectImageRegion(j,2:end)==p)==0
            for n=1:xSig
                for m=1:ySig
                    if nucleus_cluster_matrix_refresh(n,m)==p
                        NucMasks(n,m,p)=1;
                    end
                end
            end
%         end
    end
    
    seD2 = strel('disk',RingWidth);
    
    for p=1:nucleus_cluster_number_refresh
%         if sum(selectImageRegion(j,2:end)==p)==0
            DilRing = imdilate(NucMasks(:,:,p),seD2);
            RingMasks(:,:,p) = (DilRing - NucMasks(:,:,p));
            CytoArea(p)=sum(sum(RingMasks(:,:,p)));
%         end
        
    end
    
    % CYTO SIGNAL
    
    CytoSum=zeros(nucleus_cluster_number_refresh, 1);
    
    
    for p=1:nucleus_cluster_number_refresh
%         if sum(selectImageRegion(j,2:end)==p)==0
            CytoTemp=SignalImg.*RingMasks(:,:,p);
            CytoSum(p)=sum(sum(CytoTemp));
            clear CytoTemp
%         end
    end
    
    CytoSum(CytoSum==0)=[];
    CytoArea(CytoArea==0)=[];
    CytoNormArea=CytoSum./CytoArea';
    NCRatio=SigNormArea./CytoNormArea;
    
    [dataVoidRow,dataVoidCol]=size(data_out);
    dataBool=data_out==0;
    for c=1:dataVoidRow
      if dataBool(c,1)
          data_out(c,:)=[];
      end
    end
    data_out=[data_out, SigNormArea, CytoNormArea, NCRatio];
    GroupData{j}=data_out;
    close all;
end

CombinedData=GroupData{1};
for j=1:imageno-1
    CombinedData=[CombinedData; GroupData{j+1}];
end

cd(strcat(img_direc_base));
OutputDatafile=strcat('DataOutput.csv');
csvwrite(OutputDatafile, CombinedData);
close all;

clear CombinedData GroupData data_out
















