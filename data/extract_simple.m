

 startnum=0;
 runs=11;

for i = 1:runs
    runnum=startnum+i-1;

datafile = strcat('/data1/cmswank/SpinSimMirror/data/testdata.dat_',num2str(runnum));
  
% import data
    fileID = fopen(datafile);
    A = fread(fileID, 'double');
    fclose(fileID);
    Bin=A(1);
    Event=A(2);
    Var=A(3);
    A=A(4:end);
    B = reshape(A, Var, Bin, Event);
% the following are matrices (time (Bin+1), particles (Event-1))
    sx = squeeze(B(1,:,1:Event));
    sy = squeeze(B(2,:,1:Event));
    sz = squeeze(B(3,:,1:Event));
    phase=squeeze(B(4,:,1:Event));
    x =  squeeze(B(5,:,1:Event));
    y =  squeeze(B(6,:,1:Event));
    z = squeeze(B(7,:,1:Event));
    vx =  squeeze(B(8,:,1:Event));
    vy =  squeeze(B(9,:,1:Event));
    vz = squeeze(B(10,:,1:Event));
    tlarge = squeeze(B(11,:,1:Event)); 
    
    if i==1
        Sx=zeros(runs,Bin);
        Sy=zeros(runs,Bin);
        Sz=zeros(runs,Bin);
        Phase=zeros(runs,Bin);
    end
    
    if Event>1
        t = squeeze(tlarge(:,1));
        Sx(i,:)=mean(sx,2);
        Sy(i,:)=mean(sy,2);
        Sz(i,:)=mean(sz,2);
        Phase(i,:)=mean(phase,2);
   else
       t=tlarge;
       Sx=sx;
       Sy=sy;
       Sz=sz;
   end
   
   %Other stuff that might be fun to look at, or whatever, i don't know. 
   %phaseData=acos(sx1.*mean(sx1,2)+sy1.*mean(sy1,2).*sz1.*mean(sz1,2));
   %phaseNoise=std(acos(sx1.*mean(sx1,2)+sy1.*mean(sy1,2).*sz1.*mean(sz1,2))');
    %eval(['Sx',num2str(runnum),'=mean(sx1,2);']);
    %eval(['Sy',num2str(runnum),'=mean(sy1,2);']);
    %eval(['Sz',num2str(runnum),'=mean(sz1,2);']);
  
    %eval(['Signal',num2str(runnum),'=sqrt(Sx',num2str(runnum),'.^2+Sy',num2str(runnum),'.^2+Sz',num2str(runnum),'.^2);']);
end
%plotting example
plot(t,Sx(10,:));