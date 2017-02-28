# 
library(CellNOptR)
library(CNORode)
library(CNORode2017)

setwd("/Users/eduati/CPT_QSPtutorial/")

# load normalised perturbation data for LNCaP cell line in MIDAS format
Mydata<-readMIDAS(MIDASfile="perturbationData_LNCaP_MIDAS.csv")
cnolist<-makeCNOlist(Mydata, subfield=F)
cnolist$valueStimuli[cnolist$valueStimuli==0]=0.5

# load Prior Knowledge Network (PKN)
pknmodel<-readSIF("PriorKnowledgeNetwork.sif")

# compress the network to remove non-observable and non-controllable nodes
model<-preprocessing(data=cnolist, model=pknmodel, compression=TRUE, expansion=FALSE)

# set initial parameters 
ode_parameters=createLBodeContPars(model, LB_n = 1, LB_k = 0,
                                   LB_tau = 0, UB_n = 4, UB_k = 1, UB_tau = 1, default_n = 3,
                                   default_k = 0.5, default_tau = 0.01, opt_n = FALSE, opt_k = TRUE,
                                   opt_tau = TRUE, random = TRUE)

## Parameter Optimization
# essm
paramsSSm=defaultParametersSSm()
paramsSSm$local_solver = "DHC"
paramsSSm$maxtime = 36000;
paramsSSm$maxeval = Inf;
paramsSSm$atol=1e-6;
paramsSSm$reltol=1e-6;
paramsSSm$nan_fac=1000;
paramsSSm$dim_refset=30;
paramsSSm$n_diverse=1000;
paramsSSm$maxStepSize=Inf;
paramsSSm$maxNumSteps=10000;
transferFun=4;
paramsSSm$transfer_function = transferFun;

# run the optimisation algorithm
opt_pars=parEstimationLBode(cnolist,model, method="essm", ode_parameters=ode_parameters, paramsSSm=paramsSSm, lambda=0)

# plot data with optimised parameters
simulatedData=plotLBodeFitness(cnolist, model, transfer_function=transferFun, ode_parameters=opt_pars, reltol = 1e-6, atol = 1e-6, maxStepSize = 1)


