import numpy as np

num_muts_total = 16
mutations = np.array([str(x) for x in range(1,num_muts_total+1)])
mut_names = np.array(['30','35','36','57','64','65','66','79','82','83','84','85','92','95','103','113'])

num_muts_H1 = 16
H1_indices = np.arange(num_muts_total,dtype=int)
H1_mutations = mutations
H1_mut_names = mut_names

num_muts_H3 = 13
H3_required_indices = [3,8,9]
H3_required_mutations = [str(x+1) for x in H3_required_indices]
H3_remaining_indices = [0,1,2,4,5,6,7,10,11,12,13,14,15]
H3_mutations = mutations[H3_remaining_indices]
H3_mut_names = mut_names[H3_remaining_indices]

num_muts_B = 8
B_required_indices = [0,2,3,4,5,8,9,11]
B_required_mutations = [str(x+1) for x in B_required_indices]
B_remaining_indices = [1,6,7,10,12,13,14,15]
B_mutations = mutations[B_remaining_indices]
B_mut_names = mut_names[B_remaining_indices]
