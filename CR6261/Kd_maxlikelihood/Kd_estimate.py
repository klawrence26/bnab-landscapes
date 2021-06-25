import numpy as np
import scipy.optimize as opt
import scipy.stats as st

ϵ = 1e-40

def pbin(μ, σ, s, e):
    return st.norm.cdf(e, μ, σ) - st.norm.cdf(s, μ, σ)

def log_pvalue(x, μ, σ):
    res = np.zeros(x.shape)
    res[x < µ] = 2 * st.norm.cdf(x, µ, σ)[x < µ]
    res[x > µ] = 2 * (1 - st.norm.cdf(x, µ, σ)[x > µ])
    return np.mean(np.log(res + ϵ))

def gradσ_norm_pdf(x, μ, σ):
    res = st.norm.pdf(x, μ, σ)
    res[(x-μ)**2 != np.inf] *= ((x - μ)/σ)[(x-μ)**2 != np.inf]
    return res

def gradμ_pbin(μ, σ, s, e):
    return st.norm.pdf(s, μ, σ) - st.norm.pdf(e, μ, σ)

def dµµ_pbin(μ, σ, s, e):
    ress = st.norm.pdf(s, µ, σ)
    ress[(s-μ)**2 != np.inf] *= (s - μ)[(s-μ)**2 != np.inf]
    rese = st.norm.pdf(e, µ, σ)
    rese[(e-μ)**2 != np.inf] *= (e - μ)[(e-μ)**2 != np.inf]
    return (ress - rese)/(ϵ + σ**2)

def dµµ_logpbin(μ, σ, s, e):
    return (dµµ_pbin(µ, σ, s, e)*pbin(μ, σ, s, e) - gradµ_pbin(μ, σ, s, e)**2)/(ϵ + pbin(μ, σ, s, e))**2
  
    
def gradσ_pbin(μ, σ, s, e):
    return (gradσ_norm_pdf(s, μ, σ) - gradσ_norm_pdf(e, μ, σ))

def grad_logpbin(μ, σ, s, e):
    return np.row_stack([
        gradμ_pbin(μ, σ, s, e),
        gradσ_pbin(μ, σ, s, e)])/(ϵ + pbin(μ, σ, s, e))


def minimize_normal_bin(pb, lb, rb, bounds=None, µ0=None, σ0=None, ϵ=1e-40):
    """
        Given the distribution in bins pb (with bins with border lb, rb)
        Estimate both μ and σ of the normal distribution that generated them
        lb[0] = -inf
        rb[-1] = inf
        bounds: bounds[0] = (min mean, max mean), bounds[1] = (min sigma, max sigma)
        µ0, σ0: starting points on the mean and standard deviation
        Return µ, σ, fisher_information
    """

    def fun(x): # compute the likelihood
        μ, σ = x
        return -np.sum(pb * 
                      np.log(ϵ + pbin(μ, σ, lb, rb)))
    
    def grad(x):
        μ, σ = x
        return -np.sum(pb * grad_logpbin(μ, σ, lb, rb), axis=1)

    # initial values,
    bounds0 = [(lb[1] - 4, rb[-2] + 4) , (0.1, (rb[-2] - lb[1])*10)]
    if bounds is None:
        bounds = bounds0
    else:
        bounds = [(max(bounds0[0][0], bounds[0][0]), min(bounds0[0][1], bounds[0][1])),
                 (max(bounds0[1][0], bounds[1][0]), min(bounds0[1][1], bounds[1][1]))]
        
    mb = (lb + rb)/2
    if µ0 is None:
        μ0 = np.sum(mb[np.isfinite(mb)] * pb[np.isfinite(mb)])
    if σ0 is None:
        σ0 = np.max(rb[np.isfinite(rb)]) - np.min(lb[np.isfinite(lb)])
    res = opt.minimize(fun, x0=[μ0, σ0], 
                       jac=grad, 
                       bounds=bounds
                       )
    # if the bounds are reached or the minimization failed, don't trust the result
    if (res.x[0] < bounds[0][0] + 0.1) or (res.x[0] > bounds[0][1] - 0.1) or (res.success == False):
        return np.nan, np.nan, np.nan
    
    ## estimate the fisher information at that point
    fisher_info = np.sum(pb * dµµ_logpbin(res.x[0], res.x[1], lb, rb))
    return res.x[0], res.x[1], fisher_info

def log_sigmoid(x, σ, A, B):
    return np.log(A/(1 + np.exp(σ - x)) + B)


def grad_log_sigmoid(x, σ, A, B):
    x = np.array(x)
    dA = 1/(A + B + B*np.exp(σ - x))
    dB = 1/(A/(1 + np.exp(σ - x)) + B)
    dσ = -A/((1 + np.exp(σ - x)) * ((A + B)*np.exp(x - σ) + B))
    return np.column_stack([dσ, dA, dB])


def fit_log_sigmoid(x, y, sigma=None, bounds=(-np.inf, np.inf), σ0=None, A0=None, B0=None):
    """
        Find the best values for σ, A, B for the function:
        Assume A > B > 0
        \[
            log(A / (1 + exp(σ - x)) + B)
        \]
        If B is not None, keep its value fixed 
        Mean square method. 
        Return (σ, A, B)
    """
    # find reasonnable starting parameters for A & B
    # assume that we see the "full curve"
    if B0 is None:
        B0 = np.min(np.exp(y))
    if A0 is None:
        A0 = np.max(np.exp(y))
    if σ0 is None:
        σ0 = (np.min(x) + np.max(x))/2
    popt, pcov = opt.curve_fit(log_sigmoid, 
                               x,
                               y,
                               absolute_sigma=True,
                               sigma=sigma,
                               #jac=grad_log_sigmoid,
                               p0=[σ0, A0, B0],
                               bounds=bounds,
                               maxfev=400000)
    

    σ, A, B = popt[0:3]
    err = (np.sum((log_sigmoid(x, σ, A, B) - y)**2)/
                  (1e-40 + np.sum((y - y.mean())**2)))
    
    return σ, A, B, err, pcov

def sigmoid(x, σ, A, B):
    return A/(1 + np.exp(σ - x)) + B

def grad_sigmoid(x, σ, A, B):
    x = np.array(x)
    dA = 1/(1 + np.exp(σ - x))
    dB = np.ones(x.shape)
    dσ = -A*np.exp(σ - x)/(1 + np.exp(σ - x))**2
    return np.column_stack([dσ, dA, dB])

def fit_sigmoid(x, y, bounds=(-np.inf, np.inf), B=None, σ0=None):
    """
        Find the best values for σ, A, B for the function:
        Assume A > B > 0
        \[
            A / (1 + exp(σ - x)) + B
        \]
        If B is not None, keep its value fixed 
        Mean square method. 
        Return (σ, A, B)
    """
    # find reasonnable starting parameters for A & B
    # assume that we see the "full curve"
    B0 = np.min(y)
    A0 = np.max(y)
    if σ0 is None:
        σ0 = (np.min(x) + np.max(x))/2
    if B is not None:
        popt, pcov = opt.curve_fit(lambda x, σ, A: sigmoid(x, σ, A, B), 
                               x,
                               y,
                               jac=lambda x, σ, A: grad_sigmoid(x, σ, A, B)[:, :2],
                               p0=[σ0, A0],
                               bounds=[bounds[0][:2], bounds[1][:2]]
                              )    
    else:
        popt, pcov = opt.curve_fit(sigmoid, 
                               x,
                               y,
                               jac=grad_sigmoid,
                               p0=[σ0, A0, B0],
                               bounds=bounds
                              )
    σ = popt[0]
    A = popt[1]
    if B is None:
        B = popt[2]
    err = (np.sum((sigmoid(x, σ, A, B) - y)**2)/
                  np.sum((y - y.mean())**2))
    
    return σ, A, B, err


def Kd_bin_means(Rbc, bin_means, bin_stds, concentrations, cell_counts, tot_Rbc):

    B, C = Rbc.shape
    Fbc = bin_means
    pbc = Rbc * cell_counts/tot_Rbc
    pbc /= (ϵ + np.sum(pbc, axis=0))
    x = np.log(concentrations + ϵ)
    y = np.sum(pbc * Fbc, axis=0)
    yerr = np.sqrt(np.sum(pbc**2/(ϵ + np.sum(Rbc)) * Fbc**2 + pbc**2 * bin_stds**2, axis=0))
    
    okval = np.full((C,), True)
    for c in range(C):
        if np.sum(Rbc[:,c]) == 0:
            okval[c] = False
    y = y[okval]
    x = x[okval]
    yerr = yerr[okval]
    
    try:
        logKd, A, B, err, pcov =\
                fit_log_sigmoid(x, 
                                y,
                                sigma=yerr,
                                bounds=[(-50, np.min(np.exp(y))/5, np.min(np.exp(y))/5), 
                                        (-1, 2*np.max(np.exp(y)), 2*np.max(np.exp(y)))],
                                σ0=-10)
    except RuntimeError:
        return np.nan, np.nan, np.nan, np.nan, np.nan    
    return logKd, A, B, err, np.sqrt(np.diag(pcov))[0]



def Kd_gaussian_estimate(Rbc, bin_separations, concentrations,
                         cell_counts, tot_Rbc, meanlogfluo, stdlogfluo):
    """ Compute Kd.

    Arguments:
    - Rbc: numbers of reads of s in each bins, array of size B x C
    - bin_separations: fluorescence value of the bins boundaries, size B-1 x C
    - concentrations: size C, starts with 0 (no antigens)
    - cell_counts: size BxC, number of cells in each bin
    - tot_Rbc: total number of reads (every sequence) in each bins / each concentration
    Return:
    - mean of the log(fluorescence) for all concentrations
    - error on the mean of the log(fluorescence)
    - std of the log(fluorescence) for all concentrations
    - logKd, A, B
    - err: error made on the fit
    - logKderr: one standard deviation error on logKd
    """
    B, C = Rbc.shape
    Rbc = Rbc * cell_counts / tot_Rbc

    µ = np.zeros(C)
    µerr = np.zeros(C)
    σ = np.zeros(C)
    p0s = np.full(C, np.nan)
    for c in range(C):
        bin_left_edges = np.array([-np.inf] + list(np.log(bin_separations[:, c])))
        bin_right_edges = np.array(list(np.log(bin_separations[:, c])) + [np.inf])
        if np.sum(Rbc[:, c]) == 0: # no sequences
            μ[c], σ[c], μerr[c] = np.nan, np.nan, np.nan
        else: # likelihood minimisation
            pbc = Rbc[:, c]/np.sum(Rbc[:, c])
            μ[c], σ[c], fisher = minimize_normal_bin(
                pbc, bin_left_edges, bin_right_edges)
            µerr[c] = np.sqrt(
                (1/(ϵ - min(0, np.sum(Rbc[:, c])*fisher))) # fisher error
                        + 1/(ϵ + st.gmean(Rbc[:, c][Rbc[:,c]!=0]))) # error on reads
            if np.isnan(µ[c]) or np.isnan(µerr[c]):
                ## use mean-bin
                µ[c] = np.sum(pbc * meanlogfluo[:, c])
                µerr[c] = 2*np.sqrt(np.sum(pbc**2/(ϵ +  Rbc[:, c]) * meanlogfluo[:, c]**2 + pbc**2 * (stdlogfluo[:, c])**2)) # -> two standard deviation to compensate

    x = np.log(concentrations + ϵ)
    y = μ
    yerr = µerr

    ## if np.nan, remove
    ok_val = (~np.isnan(y)) & (~np.isnan(yerr))
    x = x[ok_val]
    yerr = yerr[ok_val]
    y = y[ok_val]
    
    
    if len(x) < 4:
        return µ, µerr, σ, np.nan, np.nan, np.nan, np.nan, np.nan
    
    try:
        logKd, A, B, err, pcov =\
                fit_log_sigmoid(x, 
                                y,
                                sigma=yerr,
                                bounds=[(-50, np.min(np.exp(y))/5, np.min(np.exp(y))/5), 
                                        (-1, 2*np.max(np.exp(y)), 2*np.max(np.exp(y)))],
                                σ0=-10)
    except RuntimeError:
        return µ, µerr, σ, np.nan, np.nan, np.nan, np.nan, np.nan
        
    return µ, µerr, σ, logKd, A, B, err, np.sqrt(np.diag(pcov))[0]




def Kd_estimate_bin_mean(Rbc, bin_separations, concentrations,
                         cell_counts, tot_Rbc, meanlogfluo, stdlogfluo):
    """ Compute Kd with the bin mean approach
    """
    B, C = Rbc.shape
    Rbc = Rbc * cell_counts / tot_Rbc

    µ = np.zeros(C)
    µerr = np.zeros(C)
    σ = np.zeros(C)
    p0s = np.full(C, np.nan)
    for c in range(C):
        if np.sum(Rbc[:, c]) == 0: # no sequences
            μ[c], σ[c], μerr[c] = np.nan, np.nan, np.nan
        else:
            ## use mean-bin
            pbc = Rbc[:, c]/np.sum(Rbc[:, c])
            µ[c] = np.sum(pbc * meanlogfluo[:, c])
            µerr[c] = np.sqrt(np.sum(pbc**2/(ϵ +  Rbc[:, c]) * meanlogfluo[:, c]**2 + pbc**2 * (stdlogfluo[:, c])**2))

    x = np.log(concentrations + ϵ)
    y = μ
    yerr = µerr

    ## if np.nan, remove
    ok_val = (~np.isnan(y)) & (~np.isnan(yerr))
    
    if len(x) < 4:
        return µ, µerr, yerr, np.nan, np.nan, np.nan, np.nan, np.nan
    
    
    try:
        logKd, A, B, err, pcov =\
                fit_log_sigmoid(x[ok_val], 
                                y[ok_val],
                                sigma=yerr[ok_val],
                                bounds=[(-14*np.log(10),10**3,1),(-5*np.log(10),10**5,10**5)], # angela values
                                A0=10**4,
                                B0=10,
                                σ0=-9*np.log(10),
                                # bounds=[(-50, np.min(np.exp(y))/5, np.min(np.exp(y))/5), 
                                #         (-1, 2*np.max(np.exp(y)), 2*np.max(np.exp(y)))]
                               )
    except RuntimeError:
        return µ, µerr, yerr, np.nan, np.nan, np.nan, np.nan, np.nan
        
    return µ, µerr, yerr, logKd, A, B, err, np.sqrt(np.diag(pcov))[0]