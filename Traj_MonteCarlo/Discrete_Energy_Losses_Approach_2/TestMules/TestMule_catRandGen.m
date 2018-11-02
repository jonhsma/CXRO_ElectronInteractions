%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This is a test script for the weighted categorical random generator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('..\')
nToss       =   100000;
resTest     =   cell([1 nToss]);

testCat     =   {'electron','muon','megatron','negladon','whateveron'};
testWeights =   [ 10 1 1 5 10];

for  toss   = 1:1:nToss
    resTest{toss} = weightedCategoricalRandGen(testCat,testWeights);
end

figure(44)
histogram(categorical(resTest))