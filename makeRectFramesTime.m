function frameRect=makeRectFramesTime(I,xy,z, beta, lcp)


%Description:
%This script is written for the USGS's Sandwich Beach Cam, and rectifies multiple images (I) using the specified real
%world coordinates (xy,z), the camera geometry or extrinsics (beta), and
%the camera instrinsics (lcp). The images here should already be aligned with the extrinsic calibration day. If this is not the case, look at 'auto_align.m' in order to align the imagery. The first part of the code calculates the
%shoreline elevation at the current time using buoy sea level data and runup parameterization (Stockdon, H. F., R. A. Holman, P. A. Howd, and J. Sallenger A. H. (2006)), and provides the
%rectification level (z) at the time.

%Inputs:
%1.) Add Support-Routines and UAV-Processing-Toolbox as paths. If you do
%not have these folders, you can fork them from the Coastal Imaging
%Research Networks. You also need CalcR2.m in your path which estimates
%runup. 
%2.)inputPn = file pathname that contains the images you want to rectify
%3.) Load the Sandwich buoy oceanographic data
% 4.) Load your intrinsic calibration results
% 5.) DLS_hour= make sure you converting from UTC to local time correctly
% basedon on daylight savings
% 6.) file_name= where you want to save your frameRect files 
% 7.) Hour_adj= hour adjustment based on time difference between image
% download and snap time (see note). 
% 8.) Minute_adj= minute adjustment based on time difference between image
% download and snap time (see note). 

% NOTE: This code estimates the water level based on the time difference
% between when the image was taken and when it was downloaded. These time
% differences change throughout the duration of the cameras time, and
% should be checked before processing a day's imagery.  


%Outputs:
%frameRect = rectified image structure
%frameRect.x = x coordinates [1xN]
%frameRect.y = y coordinates [1xM]
%frameRect.I = rectified image [MxNx3]

%What's next? 
%After this you can detect the shorelines using auto_shoreline_CCD.m, and
%extract the shoreline coordinates. 


%Source: Coastal Imaging Research Network GitHub: UAV-Processing-Toolbox
%Editor: Patrick Scordato, WHOI SSF 2018- United States Geological Survey-
%Woods Hole
%Date: 07/27/2018


%% Inputs. Usually only have to change the image folder inputs, and maybe the daylight savings time (convert to Greenwich Mean Time). This will change the return sea level elevation whenever converting from GMT to local time. The daylight savings adjustment is automated if you are only working on the Sandwich Beach Camera.

%1.) Add path names
addpath(genpath('D:\Scordato_SSF_2018\Source_Code\Support-Routines'))
addpath(genpath('D:\Scordato_SSF_2018\Source_Code\UAV-Processing-Toolbox'))
addpath(genpath('D:\Scordato_SSF_2018\Source_Code'))

%2.) Load images from the desired folder.
imagePn= 'D:\Scordato_SSF_2018\Projects\SandwichBeachCam\images\aligned_images\Survey_Aligned_Test_Images\March_16_2017';
d=dir(imagePn);

%3.) Oceanographic Data for the date and sea level calculations
load D:\Scordato_SSF_2018\Projects\SandwichBeachCam\oceanographic_data\Sandwidch_WL_surgeplustide_lag ;

%4.) Load in the intrinsic file
load D:\Scordato_SSF_2018\Projects\SandwichBeachCam\intrinsic_calibration\intrinsic_resullts\LCP_Results\caseLCP_720ImageSet4



%5.) Input Daylight Savings hour shifts.
% In winter, UTC= Eastern +5. In summer, UTC = Eastern +4.
%Daylight saving times in 2016 and 2017 for reference
%Summer 2016 = [2016 03 13 06 00 00];
%Winter 2016 = [2016 11 06 07 00 00];
%Summer_2017= [ 2017 03 12 06 00 00];
%Winter_2017= [2017 11 05 07 00 00];
DLS_hour= 4



%6.) Specify where you want to save your frameRect files
file_name= 'D:\Scordato_SSF_2018\Projects\SandwichBeachCam\images\frameRect_Files\All_Manual_New_RU\' ;

% 7.) Hour adjustment based on time difference between download and time
% the image was actually taken
Hour_adj= 0

% 8.) Minute adjustment based on time difference between download and time
% the image was actually taken
Minute_adj= 12

for i= 3: length(d)
    I= imread(fullfile(d(i).folder, d(i).name));
    img_name= d(i).name
 %% %Find the local time, substract download lag time, and calculate sea
 %level based on the local time
   
%     %Using image name to automatically extract estimated snap time. 
    fn=  d(i).name ;
    fn1= [erase(fn, 'L.jpg')];
    fnT= [erase(fn1, 'T')];
    space= (' ');
    fn_cell=mat2cell(fnT, 1, 2*ones(1,numel(fnT)/2));
    fn_year= strcat(fn_cell(1, 1), fn_cell(1, 2));
    fn_year2= char(fn_year);
    fn_month= char(fn_cell(1, 3));
    fn_day= char(fn_cell(1, 4));
    fn_hour= char(fn_cell(1, 5));
    fn_min= char(fn_cell(1, 6));
    fn_sec= char(fn_cell(1, 7));
    old_numTime_char= char(fn_year2, space, fn_month, space, fn_day, space, fn_hour, space, fn_min, space, fn_sec );
    
%   %User input time to get the exact time the photo was taken (manual). 
%     F= figure ;
%     imshow(I) ;
%     prompt= 'Input the time of the photo as [Year Month Day Hour Minute Seconds. Time of photo can be found in the top left corner of the image' ;
%     [old_numTime]= input(prompt) ;
%     close(F)
    
%   %Put the image time into local time and account for daylight savings
    old_numTime= transpose(str2num([old_numTime_char]));
    year= old_numTime(1, 1);
    month= old_numTime(1, 2);
    day= old_numTime(1, 3);
    hour=old_numTime(1, 4)+DLS_hour - Hour_adj;
    minute= old_numTime(1, 5) - Minute_adj;
    if minute < 0
        minute= 60+ minute ;
        hour= hour-1 ;
    end
    seconds= old_numTime(1, 6);
    numTime= [ year month day hour minute seconds ];
    
    %  Change the output file name based on the new time for saving the
    %  rectified image and image data
    img_time= num2str(numTime);
    
  
    %Calculate the water level based on the time
    tt= datenum(numTime) 
    data_num= find(T>=tt,1,'first')
    sand_total(data_num)
    wl_0= interp1(T, sand_total, tt, 'linear')
    
    %% Calculate runup 
    H = 0.2;
    T = 5;
    slope = 0.07;
    N = length(H);
    g = 9.81;
    igflag= 0 ;
    calcR2(H,T,slope,igflag);
    
    %Here is where you can change the runup value (switch between R2 and
    %R16)
    ru= ans.reduced

    %Add the runup to the sea level
    wl= ru+ wl_0
    
    %Rectification Level
    z= wl
    
    
    %%  Start the rectification process. Defines your extrinsic parameters.
   
    
    %Beta for snapshot extrinsic calibration in UTM coordinates from March 30
    %2016 Survey
    beta=  [0         0    9.0000      1.709673581931402      1.273178021597766          0]; 
    
    
    %XY Ranges and Resolution
    xy= [0 .1 400 -250 .1 100];
    
    
    %% organize indices
    %I = double(I);
    [NV,NU,NC] = size(I);
    Us = [1:NU];
    Vs = [1:NV]';
    
    %% define x,y,z grids
    x = [xy(1):xy(2): xy(3)]; y = [xy(4):xy(5): xy(6)];
    [X,Y] = meshgrid(x,y);
    
    if length(z)==1
        xyz = [X(:) Y(:) repmat(z, size(X(:)))];
    else
        xyz = [X(:) Y(:) z(:)];
    end
    
    %% Recall, Projection Matrix, P=KR[I|-C]
    %Calculate P matrix
    %define K matrix (intrinsics)
    K = [lcp.fx 0 lcp.c0U;
        0 -lcp.fy lcp.c0V;
        0  0 1];
    %define rotation matrix, R (extrinsics)
    R = angles2R(beta(4), beta(5), beta(6));
    %define identity & camera center coordinates
    IC = [eye(3) -beta(1:3)'];
    %calculate P
    P = K*R*IC ;
    %make P homogenous
    P = P/P(3,4);
    
    %% Now, convert XYZ coordinates to UV coordinates
    %convert xyz locations to uv coordinates
    UV = P*[xyz'; ones(1,size(xyz,1))];
    %homogenize UV coordinates (divide by 3 entry)
    UV = UV./repmat(UV(3,:),3,1);
    
    %convert undistorted uv coordinates to distorted coordinates
    [U,V] = distort(UV(1,:),UV(2,:),lcp);
    UV = round([U; V]);%round to the nearest pixel locations
    UV = reshape(UV,[],2); %reshape the data into something useable
    
    %find the good pixel coordinates that are actually in the image
    good = find(onScreen(UV(:,1),UV(:,2),NU,NV));
    %convert to indices
    ind = sub2ind([NV NU],UV(good,2),UV(good,1));
    
    %% Finally, grab the RGB intensities at the indices you need and fill into your XYgrid
    %preallocate final orthophoto
    Irect = zeros(length(y),length(x),3);
    
    for i = 1:NC    % cycle through R,G,B intensities
        singleBandImage = I(:,:,i); %extract the frame
        rgbIndices = singleBandImage(ind); %extract the data at the pixel locations you need
        tempImage = Irect(:,:,i); %preallocate the orthophoto size based on your x,y
        tempImage(good) = rgbIndices; %fill your extracted pixels into the frame
        Irect(:,:,i) = tempImage; %put the frame into the orthophoto
    end
    
    frameRect.x=x;
    frameRect.y=y;
    frameRect.I=Irect;
    
    %% Plot the rectified image
    rectI= figure; imagesc(frameRect.x,frameRect.y,uint8(frameRect.I));
    axis xy;axis image
    
    img_time= strcat('Date and Time: ', img_time);
    title(img_time);
    xlabel('X (m)');
    ylabel('Y (m)');
    
    
    %Save the image data (use these files to automatically extract the
    %shoreline)
    data_name= replace(img_name, '.jpg', '.mat')
    outPn= strcat(file_name, data_name)
    save(outPn, 'frameRect', 'wl');
end






