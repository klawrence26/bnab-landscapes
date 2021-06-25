import numpy as np
import csv
import sys
from sklearn.preprocessing import PolynomialFeatures
import statsmodels.api as sm

# choose statistical or biochemical epistasis
#ep_type = 'biochem' 
ep_type = 'stat'


# read in rsq data from CV folds        
num_folds = 10
if ep_type == 'stat':
	num_folds = 9
rsq_train_list_H1 = []
rsq_test_list_H1 = []
for f in range(num_folds):
	rsq_train_list_H1.append([])
	rsq_test_list_H1.append([])
	if ep_type == 'stat':
		fold = f+1
	else:
		fold = f
	with open('model_coefs/temp_9114_H1_CV_rsquared_'+str(fold)+'_'+ep_type+'.csv','r') as readfile:
		rsq_reader = csv.reader(readfile)
		header = next(rsq_reader)
		for row in rsq_reader:
			rsq_train_list_H1[f].append(float(row[1]))
			rsq_test_list_H1[f].append(float(row[2]))
        
rsq_train_list_H1 = np.array(rsq_train_list_H1)
rsq_test_list_H1 = np.array(rsq_test_list_H1)

# average over folds to get optimal order
mean_rsq_train_H1 = np.mean(rsq_train_list_H1,axis=0)
stdev_rsq_train_H1 = np.std(rsq_train_list_H1,axis=0)
mean_rsq_test_H1 = np.mean(rsq_test_list_H1,axis=0)
stdev_rsq_test_H1 = np.std(rsq_test_list_H1,axis=0)

optimal_H1_order = np.argmax(mean_rsq_test_H1)
print('Optimal order, H1: ',optimal_H1_order,file=sys.stdout,flush=True)

print(mean_rsq_test_H1)

# write averages to new file
with open('model_coefs/9114_H1_CV_rsquared_'+ep_type+'.csv','w') as writefile:
	rsq_writer = csv.writer(writefile)
	rsq_writer.writerow(['Optimal order: '+str(optimal_H1_order)])
	rsq_writer.writerow(['Type','Order','Mean','Std'])
	for i in range(len(mean_rsq_train_H1)):
		rsq_writer.writerow(['Train',str(i),mean_rsq_train_H1[i],stdev_rsq_train_H1[i]])
	for i in range(len(mean_rsq_test_H1)):
		rsq_writer.writerow(['Test',str(i),mean_rsq_test_H1[i],stdev_rsq_test_H1[i]])
	writefile.close()



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
    

# # Fit final models

np.random.seed(2112)
indices_permuted_H1 = np.random.permutation(np.arange(len(genos_H1)))

# fit models of increasing order
for order in range(1,optimal_H1_order+1):

    genos_H1_permuted = genos_H1[indices_permuted_H1]
    phenos_H1_permuted = phenos_H1[indices_permuted_H1]
    print('Order: ',str(order),file=sys.stdout,flush=True)
    poly_H1_current = PolynomialFeatures(order,interaction_only=True)
    genos_H1_current = poly_H1_current.fit_transform(genos_H1_permuted)

    # fit
    reg_H1_current = sm.OLS(phenos_H1_permuted,genos_H1_current).fit()
    reg_H1_coefs_current = reg_H1_current.params
    reg_H1_CIs_current = reg_H1_current.conf_int(alpha=0.05/float(len(reg_H1_coefs_current)), cols=None)
    reg_H1_stderr = reg_H1_current.bse
    reg_H1_pvalues = reg_H1_current.pvalues
    
    num_sig = len(np.where(reg_H1_pvalues < 0.05/float(len(reg_H1_coefs_current)))[0])

    predicted_phenos_permuted_H1 = reg_H1_current.predict(genos_H1_current)
    rsquared_H1_current = reg_H1_current.rsquared
    print('Params: ',len(reg_H1_coefs_current),file=sys.stdout,flush=True)
    print('Performance: ',rsquared_H1_current,file=sys.stdout,flush=True)
    print(num_sig,file=sys.stdout,flush=True)
	 
    #sys.exit()

    # write model to file
    coef_names = poly_H1_current.get_feature_names(input_features = mutations_H1)
    with open('model_coefs/9114_H1_'+str(order)+'order_'+ep_type+'.txt','w') as writefile:
        coef_writer = csv.writer(writefile,delimiter='\t')
        coef_writer.writerow(['Params: ',len(reg_H1_coefs_current)])
        coef_writer.writerow(['Performance: ',rsquared_H1_current])
        coef_writer.writerow(['Term','Coefficient','Standard Error','p-value','95% CI lower','95% CI upper'])
        coef_writer.writerow(['Intercept',reg_H1_coefs_current[0]])
        for i in range(1,len(reg_H1_coefs_current)):
            coef_writer.writerow([','.join(coef_names[i].split(' ')),reg_H1_coefs_current[i],reg_H1_stderr[i],
                                  reg_H1_pvalues[i],reg_H1_CIs_current[i][0],reg_H1_CIs_current[i][1]])
        writefile.close()


