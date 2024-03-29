clc
clear

  data1=readmatrix("Trial _trial_ID1_Ch2_Channel 2.csv");
  data2=readmatrix("Trial _trial_ID1_Ch3_Channel 3.csv");
  data3=readmatrix("Trial _trial_ID1_Ch5_Channel 5.csv");   
  [Hm01,f1,S1,zsub1,Tp1]=spectrumWaves(data1(:,3),data2(:,3),data2(:,3),0.01);

  