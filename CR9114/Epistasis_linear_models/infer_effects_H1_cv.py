import numpy as np
import csv
import sys
from sklearn.preprocessing import PolynomialFeatures
import statsmodels.api as sm


# choose statistical or biochemical epistasis
#ep_type = 'biochem' 
ep_type = 'stat'



# read in data

geno_vectors_H1 = []
phenos_H1 = []

mutations_H1 = [str(x) for x in range(1,17)]


with open('../Kd_meanbin/kd_processed/20210323_9114_HA_unadj_fil_merg.csv','r') as readfile:
    kd_reader = csv.reader(readfile)
    header = next(kd_reader)
    for row in kd_reader:
        geno = row[0]
        
        geno_vec = np.array([float(x) for x in geno])

        pheno_H1 = row[4]
        
        if len(pheno_H1) != 0:  
            geno_vectors_H1.append(geno_vec)
            phenos_H1.append(float(pheno_H1))
    readfile.close()

phenos_H1 = np.array(phenos_H1)

genos_H1 = np.empty((len(phenos_H1),len(geno_vectors_H1[0])))
for i in range(len(phenos_H1)):
    genos_H1[i] = geno_vectors_H1[i][:]
    
if ep_type == 'stat':
    genos_H1 = 2*(genos_H1-0.5)

print(genos_H1.shape,phenos_H1.shape)


# # CV to choose optimal order

f = int(sys.argv[1])
num_folds = 1
max_order = 6

# set up permutation
np.random.seed(2112)
indices_permuted_H1 = np.random.permutation(np.arange(len(genos_H1)))
size_test_H1 = int(0.1*len(genos_H1))
size_train_H1 = len(genos_H1)-size_test_H1
print(size_test_H1,size_train_H1,file=sys.stdout,flush=True)

# lists to store r squared values
rsq_train_list_H1 = np.zeros(max_order+1)
rsq_test_list_H1 = np.zeros(max_order+1)



# loop over CV folds
for fold in range(1):

    # get train & test sets
    start = int(f*size_test_H1)
    stop = int((f+1)*size_test_H1)
    genos_train_H1 = np.concatenate((genos_H1[indices_permuted_H1[:start]],genos_H1[indices_permuted_H1[stop:]]))
    genos_test_H1 = genos_H1[indices_permuted_H1[start:stop]]
    phenos_train_H1 = np.concatenate((phenos_H1[indices_permuted_H1[:start]],phenos_H1[indices_permuted_H1[stop:]]))
    phenos_test_H1 = phenos_H1[indices_permuted_H1[start:stop]]
    
    print('Fold: ',f,file=sys.stdout,flush=True)
    
    # initialize zero-order (intercept-only) model
    print('Order: 0',file=sys.stdout,flush=True)
    genos_train_H1_previous = np.full(len(genos_train_H1),1.0)
    genos_test_H1_previous = np.full(len(genos_test_H1),1.0)

    reg_H1_previous = sm.OLS(phenos_train_H1,genos_train_H1_previous).fit()
    reg_H1_coefs_previous = reg_H1_previous.params

    rsquared_train_H1_previous = reg_H1_previous.rsquared
    rsquared_test_H1_previous = 1-np.sum((phenos_test_H1-reg_H1_previous.predict(genos_test_H1_previous))**2)/np.sum((phenos_test_H1-np.mean(phenos_test_H1))**2)
    rsq_train_list_H1[0] = rsquared_train_H1_previous
    rsq_test_list_H1[0] = rsquared_test_H1_previous

    mean_pheno_train = np.mean(phenos_train_H1)
    mean_pheno_test = np.mean(phenos_test_H1)


    # fit models of increasing order
    for order in range(1,max_order+1):
        print('Order: ',str(order),file=sys.stdout,flush=True)
        poly_H1_current = PolynomialFeatures(order,interaction_only=True)
        genos_train_H1_current = poly_H1_current.fit_transform(genos_train_H1)
        genos_test_H1_current = poly_H1_current.fit_transform(genos_test_H1)

        reg_H1_current = sm.OLS(phenos_train_H1, genos_train_H1_current).fit()
        reg_H1_coefs_current = reg_H1_current.params
        reg_H1_CIs_current = reg_H1_current.conf_int(alpha=0.05, cols=None)
        reg_H1_stderr = reg_H1_current.bse
    
        rsquared_train_H1_current = reg_H1_current.rsquared
        rsquared_test_H1_current = 1-np.sum((phenos_test_H1-reg_H1_current.predict(genos_test_H1_current))**2)/np.sum((phenos_test_H1-np.mean(phenos_test_H1))**2)
        rsq_train_list_H1[order] = rsquared_train_H1_current
        rsq_test_list_H1[order] = rsquared_test_H1_current
                
        
print(rsq_train_list_H1,file=sys.stdout,flush=True)
print(rsq_test_list_H1,file=sys.stdout,flush=True)


# write r2 values to file
with open('model_coefs/temp_9114_H1_CV_rsquared_'+str(f)+'_'+ep_type+'.csv','w') as writefile:
	rsq_writer = csv.writer(writefile)
	rsq_writer.writerow(['Order','Train R2','Test R2'])
	for i in range(len(rsq_train_list_H1)):
		rsq_writer.writerow([i,rsq_train_list_H1[i],rsq_test_list_H1[i]])
	writefile.close()

