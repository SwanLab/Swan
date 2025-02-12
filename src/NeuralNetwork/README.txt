This folder contains all the matlab functions/classes created to develope the thesis.
It also contains the "datasets.mat" which is a vector of strings containing the name 
of the datasets for convenience. 

The initialisation of the 3 main classes have to be summoned as follows:

***For class Data***
data = Data(fn,tr,p)

where fn is the name of the file where is the datasets (.csv), tr is the ratio used for test
(30 means 30%), and p is the degree of polynomial augmentation to the data desired.

***For class Network***
It is possible to give 2 or 6 input arguments

network = Network(data,structure)
network = Network(data,structure,costfunction,hiddenunit,outputunit,lambda)

Where data is an object Data with the data wanted to train the ANN with, structure is a vector
containing the number of neurons at each layer (input and output layer also), costfunction is 
the cost function to use (options: '-loglikelihood' or 'L2'), hiddenunit is the activation 
function desired in the hidden layers (options: 'sigmoid', 'tanh' and 'ReLU'), outputunit is 
the activation function desired in the output layer (options: 'sigmoid' and 'softmax') and
lambda is a real number (value of the L2 regularization parameter).

With 2 inputs the default parameters are:
network = Network(data,structure,'-loglikelihood','ReLU','softmax',0);

***For class Trainer***
It is possible to give 4, 7 or 8 input arguments

optimizer = Trainer.create(network,algorithm,learningRate,momentum);
optimizer = Trainer.create(network,algorithm,learningRate,momentum,batch,options,linesearch);
optimizer = Trainer.create(network,algorithm,learningRate,momentum,batch,options,linesearch,nplt);

Where network is the network which has to be trained, algorithm is the algorithm used (options: 
'SGD', 'Nesterov', 'RMSProp'), learningRate is a real number, momentum is a real number (alpha
for Nesterov and rho for RMSProp), batch is an integer (batch size), options are a set of minimization
criteria, linesearch is the type of linesearch (options: 'static', 'decay', 'dynamic', 'fminbnd') and 
nplt is a real number if the plot of the minimization is desired which represent every how many epochs
to plot.

With 4 inputs the default parameters are:
optimizer = Trainer.create(network,algorithm,learningRate,momentum,200,opt,'static',0);
opt.optTolerance  = 1*10^-6;    opt.maxevals      = 1000;
opt.maxepochs     = 1000   ;    opt.earlyStop     = 1000;
opt.time          = Inf([1,1]); opt.fv         = 10^-4;
