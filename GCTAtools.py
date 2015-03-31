## Module of python functions to help with the processing and analysis of GCTA results files
## Davis McCarthy, July 2014

## Import a whole bunch of packages and modules used
import re
import pandas as pd
import numpy as np
import scipy.stats as st
from struct import unpack, calcsize


#def parseHE(hefile, names_var_comps = ['Enhancer', 'Other'], npcs = 10, minMAF = '0.001', ldpruning = 'No pruning', 
#             relthresh = 'relmax0.05', cohort = 'ALL', grmalg = 's = -1', fittype = 'HE', constraint = 'Unconstrained', 
#             remlalg = 'NA'):
#    """ 
#    returns the variance components results and parameter covariance matrix from a GEAR *.he file 
#    Davis McCarthy, October 2014
#    """
#    with open ( hefile , 'r' ) as IN:
#        results = pd.DataFrame(columns=('Coef', 'Estimate', 'se'))
#        variance_components = pd.DataFrame(columns=('Component', 'Variance', 'SE'))
#        variance_explained_liability = pd.DataFrame(columns=('Source', 'Variance', 'SE'))
#        variance_matrix = []
#        varcovar = False
#        results_found = False
#        sum_ve_est = None
#        sum_ve_se = None
#        for line in IN:
#            line = line.strip()
#            ## If line starts with 'b' then take it as one of the coefficients
#            if re.match('^b', line):
#                results_found = True
#                line = line.strip()
#                line_list = re.split('\s', line)
#                if( len(line_list) > 1 ):
#                    row = {'Coef': line_list[0], 'Estimate': line_list[1], 'se': line_list[2]}
#                    results = results.append(row, ignore_index=True)
#            elif( re.match('^variance-covariance', line) ):
#                varcovar = True
#                continue
#            if( varcovar and re.match('$', line) ):
#                varcovar = False
#                continue
#            if( varcovar ):
#                line_list = re.split('\s', line)
#                variance_matrix.append(line_list)
#    ## If logfile is complete and results are found, then process and format the results
#    if( results_found ):
#        ## Sort the covariance matrix for variance components and compute enrichment scores
#        n = len(variance_matrix)
#        mtrx = [[0 for x in range(n)] for x in range(n)]
#        for i in range(n):
#            for j in range(n):
#                if( i >= j ): 
#                    mtrx[i][j] = variance_matrix[i][j]
#                else: 
#                    mtrx[i][j] = variance_matrix[j][i]
#        covar_matrix = np.array(mtrx, np.float64)
#        var_comps = np.array(variance_components['Variance'], np.float64)
#        enrichment = calcEnrichmentScore(var_comps, covar_matrix)
#        ## Manage metadata
#        ncomps = len(names_var_comps)
#        indices = variance_explained_liability.index 
#        if( len(minMAF) == 1 ):
#            minMAF = [minMAF for i in range(ncomps)]
#        metadata = pd.DataFrame({'minMAF': pd.Series(minMAF, index=indices),
#                                 'Var_Comp': pd.Series(names_var_comps, index=indices),
#                                 'Sum_VG_Est': pd.Series([sum_ve_est for i in range(ncomps)], index=indices),
#                                 'Sum_VG_SE': pd.Series([sum_ve_se for i in range(ncomps)], index=indices),
#                                 'LD_Pruning': pd.Series([ldpruning for i in range(ncomps)], index=indices),
#                                 'N_PCs': pd.Series([npcs for i in range(ncomps)], index=indices),
#                                 'Rel_Thresh': pd.Series([relthresh for i in range(ncomps)], index=indices),
#                                 'Cohort': pd.Series([cohort for i in range(ncomps)], index=indices),
#                                 'GRM_Alg': pd.Series([grmalg for i in range(ncomps)], index=indices),
#                                 'REML_Alg': pd.Series([remlalg for i in range(ncomps)], index=indices),
#                                 'Fit_Type': pd.Series([fittype for i in range(ncomps)], index=indices),
#                                 'Constraint': pd.Series([constraint for i in range(ncomps)], index=indices)
#                                 })
#        enrich_df = enrichment['enrichment']
#        enrich_df.index = indices
#        enrich_df.columns = ['Enrichment', 'Enrich_SE']
#        ve_out = variance_explained_liability.join(enrich_df, how='inner')
#        ve_out = ve_out.join(metadata, how='inner')
#        options = None
#        out = {'enrichment': enrichment, 'options': options, 'results': results, 'variance_components': variance_components, 
#           'variance_explained': ve_out, 'covar_matrix': covar_matrix}
#    else:
#        out = {'enrichment': None, 'options': options, 'results': None, 'variance_components': None, 
#           'variance_explained': None, 'covar_matrix': None}
#    return out
#    
#    
#def calcVEfromHE(coefs, covar_matrix):
#    """
#    returns VE/heritability estimates from the coefficients from HE regression
#    """
#    n = len(coefs)
#    ## Compute h2 ests
#    genet_coefs = coefs[1:]
#    b0 = coefs[0]
#    h2_ests = [-1 * b / b0 for b in genet_coefs]
#    ## Compute grad matrix for component heritability estimates
#    grad_h2 = [[0 for x in range(n-1)] for x in range(n)]
#    for i in range(n-1):
#        for j in range(n):
#            if( j == 0 ):
#                k = i + 1
#                grad_h2[i][j] = -coefs[k] * np.log(b0)
#            elif( i == j ): 
#                grad_h2[i][j] = -1 / b0
#            else: 
#                grad_h2[i][j] = 0
#    grad_h2 = np.matrix(grad_h2, np.float64)
#    ## Compute covariance matrix for h2 ests by delta method
#    h2_var = grad_h2 * covar_matrix * grad_h2.transpose()
#    ## Get standard errors for h2 ests
#    h2_se_ests = np.sqrt(np.diag(h2_var))    
#    out = {'h2_ests': h2_ests, 'se': h2_se_ests}
#    return out
#    
#    
#def calcVETotalfromHE(coefs, covar_matrix):
#    """
#    returns total heritability estimate from a set of HE coefficients and 
#    covariance matrix
#    """
#    n = len(coefs)
#    b0 = coefs[0]
#    bi_total = sum(coefs[1:])
#    h2_total = -1 * bi_total / b0
#    ## Compute se's for total (sum) h2 est
#    grad_h2_total = np.append([-np.log(b0) * bi_total], [(-1 / b0) for x in range(n-1)])
#    grad_h2_total = np.matrix(grad_h2_total)
#    h2_total_var = grad_h2_total * covar_matrix * grad_h2_total.transpose()
#    h2_total_se_est = np.sqrt(h2_total_var)
#    out = {'h2_total': h2_total, 'se': h2_total_se_est}
#    return out


def parseLog(logfile, names_var_comps = ['Enhancer', 'Other'], npcs = 10, minMAF = '0.001', ldpruning = 'No pruning', 
             relthresh = 'relmax0.05', cohort = 'IP_ALL', grmalg = 's = -1', fittype = 'Joint', constraint = 'Constrained', 
             remlalg = 'EM'):
    """ 
    returns the GCTA options, variance components results and parameter covariance matrix from a GCTA log file 
    Davis McCarthy, June 2014
    """
    with open ( logfile , 'r' ) as IN:
        options = pd.DataFrame(columns=('Option', 'Arg_Supplied'))
        results = pd.DataFrame(columns=('Source', 'Variance', 'SE'))
        variance_components = pd.DataFrame(columns=('Source', 'Variance', 'SE'))
        variance_explained_liability = pd.DataFrame(columns=('Source', 'Variance', 'SE'))
        variance_matrix = []
        varcovar = False
        results_found = False
        sum_ve_est = None
        sum_ve_se = None
        for line in IN:
            line = line.strip()
            ## If line starts with '--' then take it as the list of options to GCTA and save to dict
            if re.match('^--', line):
                line = line.strip('-')
                line_list = re.split('\s', line)
                if( len(line_list) > 1 ):
                    row = {'Option': line_list[0], 'Arg_Supplied': line_list[1]}
                    options = options.append(row, ignore_index=True)
            elif( re.match('V[()Gep0-9]+', line) or re.match('Sum of V[()Gep0-9]+', line) or re.match('Vp', line) ):
                results_found = True
                line_list = re.split('\t', line)
                row = {'Source': line_list[0], 'Variance': line_list[1], 'SE': line_list[2]}
                results = results.append(row, ignore_index=True)
                if( re.match('V[()Ge0-9]+\s', line) ):
                    variance_components = variance_components.append(row, ignore_index=True)
                if( re.match('V[()Ge0-9]+/Vp_L', line) ):
                    variance_explained_liability = variance_explained_liability.append(row, ignore_index=True)
                if( re.match('Sum of V[(G)]+_L', line) ):
                    sum_ve_est = row['Variance']
                    sum_ve_se = row['SE']
            if( re.match('^Variance/Covariance Matrix', line) ):
                varcovar = True
                continue
            if( varcovar and re.match('$', line) ):
                varcovar = False
                continue
            if( varcovar ):
                line_list = re.split('\s', line)
                variance_matrix.append(line_list)
    ## If logfile is complete and results are found, then process and format the results
    if( results_found ):
        ## Sort the covariance matrix for variance components and compute enrichment scores
        n = len(variance_matrix)
        mtrx = [[0 for x in range(n)] for x in range(n)]
        for i in range(n):
            for j in range(n):
                if( i >= j ): 
                    mtrx[i][j] = variance_matrix[i][j]
                else: 
                    mtrx[i][j] = variance_matrix[j][i]
        covar_matrix = np.array(mtrx, np.float64)
        var_comps = np.array(variance_components['Variance'], np.float64)
        enrichment = calcEnrichmentScore(var_comps, covar_matrix)
        ## Manage metadata
        ncomps = len(names_var_comps)
        indices = variance_explained_liability.index 
        if( len(minMAF) == 1 ):
            minMAF = [minMAF for i in range(ncomps)]
        metadata = pd.DataFrame({'minMAF': pd.Series(minMAF, index=indices),
                                 'Var_Comp': pd.Series(names_var_comps, index=indices),
                                 'Sum_VG_Est': pd.Series([sum_ve_est for i in range(ncomps)], index=indices),
                                 'Sum_VG_SE': pd.Series([sum_ve_se for i in range(ncomps)], index=indices),
                                 'LD_Pruning': pd.Series([ldpruning for i in range(ncomps)], index=indices),
                                 'N_PCs': pd.Series([npcs for i in range(ncomps)], index=indices),
                                 'Rel_Thresh': pd.Series([relthresh for i in range(ncomps)], index=indices),
                                 'Cohort': pd.Series([cohort for i in range(ncomps)], index=indices),
                                 'GRM_Alg': pd.Series([grmalg for i in range(ncomps)], index=indices),
                                 'REML_Alg': pd.Series([remlalg for i in range(ncomps)], index=indices),
                                 'Fit_Type': pd.Series([fittype for i in range(ncomps)], index=indices),
                                 'Constraint': pd.Series([constraint for i in range(ncomps)], index=indices)
                                 })
        enrich_df = enrichment['enrichment']
        enrich_df.index = indices
        enrich_df.columns = ['Enrichment', 'Enrich_SE']
        ve_out = variance_explained_liability.join(enrich_df, how='inner')
        ve_out = ve_out.join(metadata, how='inner')
        out = {'enrichment': enrichment, 'options': options, 'results': results, 'variance_components': variance_components, 
           'variance_explained': ve_out, 'covar_matrix': covar_matrix}
    else:
        out = {'enrichment': None, 'options': options, 'results': None, 'variance_components': None, 
           'variance_explained': None, 'covar_matrix': None}
    return out



def parseLogContTrait(logfile, names_var_comps = ['Enhancer', 'Other'], npcs = 10, minMAF = '0.001', ldpruning = 'No pruning', 
             relthresh = 'relmax0.05', cohort = 'IP_ALL', grmalg = 's = -1', fittype = 'Joint', constraint = 'Constrained', 
             remlalg = 'EM'):
    """ 
    returns the GCTA options, variance components results and parameter covariance matrix from a GCTA log file 
    Davis McCarthy, November 2014
    """
    with open ( logfile , 'r' ) as IN:
        options = pd.DataFrame(columns=('Option', 'Arg_Supplied'))
        results = pd.DataFrame(columns=('Source', 'Variance', 'SE'))
        variance_components = pd.DataFrame(columns=('Source', 'Variance', 'SE'))
        variance_explained = pd.DataFrame(columns=('Source', 'Variance', 'SE'))
        variance_matrix = []
        varcovar = False
        results_found = False
        sum_ve_est = None
        sum_ve_se = None
        for line in IN:
            line = line.strip()
            ## If line starts with '--' then take it as the list of options to GCTA and save to dict
            if re.match('^--', line):
                line = line.strip('-')
                line_list = re.split('\s', line)
                if( len(line_list) > 1 ):
                    row = {'Option': line_list[0], 'Arg_Supplied': line_list[1]}
                    options = options.append(row, ignore_index=True)
            elif( re.match('V[()Gep0-9]+', line) or re.match('Sum\s+of\s+V', line) or re.match('Vp', line) ):
                results_found = True
                line_list = re.split('\t', line)
                row = {'Source': line_list[0], 'Variance': line_list[1], 'SE': line_list[2]}
                results = results.append(row, ignore_index=True)
                if( re.match('V[()Ge0-9]+\s', line) ):
                    variance_components = variance_components.append(row, ignore_index=True)
                if( re.match('V[()Ge0-9]+/Vp\s', line) ):
                    variance_explained = variance_explained.append(row, ignore_index=True)
                if( re.match('Sum\s+of', line) ):
                    sum_ve_est = row['Variance']
                    sum_ve_se = row['SE']
            if( re.match('^Variance/Covariance Matrix', line) ):
                varcovar = True
                continue
            if( varcovar and re.match('$', line) ):
                varcovar = False
                continue
            if( varcovar ):
                line_list = re.split('\s', line)
                variance_matrix.append(line_list)
    ## If logfile is complete and results are found, then process and format the results
    if( results_found ):
        ## Sort the covariance matrix for variance components and compute enrichment scores
        n = len(variance_matrix)
        mtrx = [[0 for x in range(n)] for x in range(n)]
        for i in range(n):
            for j in range(n):
                if( i >= j ): 
                    mtrx[i][j] = variance_matrix[i][j]
                else: 
                    mtrx[i][j] = variance_matrix[j][i]
        covar_matrix = np.array(mtrx, np.float64)
        var_comps = np.array(variance_components['Variance'], np.float64)
        enrichment = calcEnrichmentScore(var_comps, covar_matrix)
        ## Manage metadata
        ncomps = len(names_var_comps)
        indices = variance_explained.index 
        if( len(minMAF) == 1 ):
            minMAF = [minMAF for i in range(ncomps)]
        metadata = pd.DataFrame({'minMAF': pd.Series(minMAF, index=indices),
                                 'Var_Comp': pd.Series(names_var_comps, index=indices),
                                 'Sum_VG_Est': pd.Series([sum_ve_est for i in range(ncomps)], index=indices),
                                 'Sum_VG_SE': pd.Series([sum_ve_se for i in range(ncomps)], index=indices),
                                 'LD_Pruning': pd.Series([ldpruning for i in range(ncomps)], index=indices),
                                 'N_PCs': pd.Series([npcs for i in range(ncomps)], index=indices),
                                 'Rel_Thresh': pd.Series([relthresh for i in range(ncomps)], index=indices),
                                 'Cohort': pd.Series([cohort for i in range(ncomps)], index=indices),
                                 'GRM_Alg': pd.Series([grmalg for i in range(ncomps)], index=indices),
                                 'REML_Alg': pd.Series([remlalg for i in range(ncomps)], index=indices),
                                 'Fit_Type': pd.Series([fittype for i in range(ncomps)], index=indices),
                                 'Constraint': pd.Series([constraint for i in range(ncomps)], index=indices)
                                 })
        enrich_df = enrichment['enrichment']
        enrich_df.index = indices
        enrich_df.columns = ['Enrichment', 'Enrich_SE']
        ve_out = variance_explained.join(enrich_df, how='inner')
        ve_out = ve_out.join(metadata, how='inner')
        out = {'enrichment': enrichment, 'options': options, 'results': results, 'variance_components': variance_components, 
           'variance_explained': ve_out, 'covar_matrix': covar_matrix}
    else:
        out = {'enrichment': None, 'options': options, 'results': None, 'variance_components': None, 
           'variance_explained': None, 'covar_matrix': None}
    return out

    

def parseHsq(hsqfile, names_var_comps = ['Enhancer', 'Other'], npcs = 10, minMAF = 0.001, ldpruning = 'No pruning', 
             relthresh = 'relmax0.05', cohort = 'IP_ALL', grmalg = 's = -1', fittype = 'Joint', constraint = 'Constrained'):
    """
    return relevant results and fit details from an Hsq file produced by GCTA
    Davis McCarthy, June 2014
    """
    with open ( hsqfile , 'r' ) as IN:
        fit = pd.DataFrame(columns=('Attribute', 'Value'))
        results = pd.DataFrame(columns=('Source', 'Variance', 'SE'))
        variance_explained_liability = pd.DataFrame(columns=('Source', 'Variance', 'SE'))
        for line in IN:
            line = line.strip()
            if( re.match('V[()Gep0-9]+', line) or re.match('Sum of V[()Gep0-9]+', line) or re.match('Vp', line) ):
                line_list = re.split('\t', line)
                row = {'Source': line_list[0], 'Variance': line_list[1], 'SE': line_list[2]}
                results = results.append(row, ignore_index=True)
                if( re.match('V[()Ge0-9]+/Vp_L', line) ):
                    variance_explained_liability = variance_explained_liability.append(row, ignore_index=True)
            elif( re.match('[lLdPn]+', line) ):
                line_list = re.split('\t', line)
                row = {'Attribute': line_list[0], 'Value': line_list[1]}
                fit = fit.append(row, ignore_index=True)
    ncomps = len(names_var_comps)
    indices = variance_explained_liability.index 
    metadata = pd.DataFrame({'minMAF': pd.Series([minMAF for i in range(ncomps)], index=indices),
                             'Var_Comp': pd.Series(names_var_comps, index=indices),
                             'LD_Pruning': pd.Series([ldpruning for i in range(ncomps)], index=indices),
                             'N_PCs': pd.Series([npcs for i in range(ncomps)], index=indices),
                             'Rel_Thresh': pd.Series([relthresh for i in range(ncomps)], index=indices),
                             'Cohort': pd.Series([cohort for i in range(ncomps)], index=indices),
                             'GRM_Alg': pd.Series([grmalg for i in range(ncomps)], index=indices),
                             'Fit_Type': pd.Series([fittype for i in range(ncomps)], index=indices),
                             'Constraint': pd.Series([constraint for i in range(ncomps)], index=indices)
                             })
    fit_expand = fit.T[1:]
    fit_expand.columns = fit['Attribute']
    fit_expand.drop(['LRT', 'Pval'], inplace = True, axis = 1)
    fitrow = fit_expand.iloc[0]
    for i in range(0, (ncomps - 1)):
        fit_expand = fit_expand.append(fitrow, ignore_index = True)
    fit_expand.index = indices
    ve_out = variance_explained_liability.join(metadata, how='inner')
    ve_out = ve_out.join(fit_expand, how='inner')
    out = {'results': results, 'variance_explained': ve_out}
    return out
    
def addEnrichmentResults(df, num_vars, sum_num_vars = None):
    """
    convenience function to add columns with enrichment results to a dataframe produced by parseLog()
    nvars should be an array of integers/floats
    Davis McCarthy, June 2014
    """
    if sum_num_vars is None:
        sum_num_vars = sum(num_vars)
    df['N_Vars'] = num_vars  
    df['Prop_Vars'] = num_vars / sum_num_vars
    diff_enrich_expected = (df['Enrichment'] - df['Prop_Vars'])
    fold_enrichment =  df['Enrichment'] / df['Prop_Vars']
    fold_enrichment[ df['Enrichment'] < 0 ] = -fold_enrichment[ df['Enrichment'] < 0 ]
    df['Fold_Enrichment'] = fold_enrichment
    df['Enrich_ZScore'] = (diff_enrich_expected / df['Enrich_SE'])
    df['Enrich_Pval'] = 2*st.norm.sf(abs(df['Enrich_ZScore']))
    cols = ['Var_Comp', 'Source', 'Variance', 'SE', 'Sum_VG_Est', 'Sum_VG_SE', 'Enrichment', 'Fold_Enrichment', 'Enrich_SE', 'Enrich_ZScore', 'Enrich_Pval', 'minMAF',
            'N_Vars', 'Prop_Vars', 'Cohort', 'Constraint', 'Fit_Type', 'GRM_Alg', 'REML_Alg', 'LD_Pruning', 'N_PCs', 'Rel_Thresh']
    return df[cols]


def calcEnrichmentScore(var_comps, covar_matrix):
    """
    given variance component estimates and a covariance matrix it computes the enrichment score (percentage of genetic
    variance explained by a component)
    assumes that the final variance component is the error variance component, which is not used (but expected for completeness)
    returns a dataframe with enrichment scores and standard errors computed using the delta method
    Davis McCarthy, June 2014
    """
    n = len(var_comps)
    ## Compute enrichment score and grad fn
    enrich_scores = calcGradofEnrichmentFn(var_comps)
    ## Pass grad matrix and covariance matrix (of only the genetic variance components)
    variance_enrich_scores = calcVarianceEnrichScores(enrich_scores['grad_scores'], covar_matrix[:(n-1), :(n-1)])
    ## Compute SEs for the enrichment scores
    se = np.sqrt(np.diag(variance_enrich_scores))
    ## Organise results for output
    enrich_df = pd.DataFrame({'Enrichment_Score': enrich_scores['enrichment_scores'], 'SE': se})
    out = {'enrichment': enrich_df, 'variance_enrichment': variance_enrich_scores, 'grad_enrichment': enrich_scores['grad_scores']}
    return out

    
def calcGradofEnrichmentFn(var_comps):
    """
    return the enrichment scores as a vector (1-d array) and the grad of the enrichment function (n-1 x n-1 matrix)
    Davis McCarthy, June 2014
    """
    n = len(var_comps)
    ## Compute enrichment scores
    genet_var_comps = var_comps[:(n-1)]
    sum_var_comps = sum(genet_var_comps)
    sq_sum_var_comps = sum_var_comps ** 2
    scores = [x / sum_var_comps for x in genet_var_comps]
    ## Compute grad matrix of enrichment scores
    grad_scores = [[0 for x in range(n-1)] for x in range(n-1)]
    for i in range(n-1):
        for j in range(n-1):
            if( i == j ): 
                grad_scores[i][j] = sum(np.delete(genet_var_comps, i)) / sq_sum_var_comps
            else: 
                grad_scores[i][j] = -genet_var_comps[i] / sq_sum_var_comps
    grad_scores = np.matrix(grad_scores, np.float64)
    out = {'enrichment_scores': scores, 'grad_scores': grad_scores}
    return out
    
    
def calcVarianceEnrichScores(grad_scores, covar_matrix):
    """
    return SEs for enrichment scores
    Davis McCarthy, June 2014
    """
    ## Compute standard errors for enrichment scores
    return grad_scores * covar_matrix * grad_scores.transpose()
    

def showEnrichment(parsed_log):
    print(parsed_log['enrichment']['enrichment'])
    
    
def ReadGRMBinN(prefix):
    """
    read the number of SNPs used to compute a GRM from GCTA-format binary file
    adapted from an R function on GCTA website
    Davis McCarthy, June 2014
    """
    NFileName = prefix + ".grm.N.bin"
    entry_format = 'f' #N is stored as a float in the binary file
    entry_size = calcsize(entry_format)
    ## Open the binary file
    with open(NFileName, mode='rb') as f:
        #entry_count = os.fstat(f.fileno()).st_size / entry_size
        record = f.read(entry_size)
        N = unpack(entry_format, record)[0]
        N = int(N)
    return(N)


def sum_n_vec(n):
    """
    return a vector giving the one less than the sum of the first i
    integers for i from 1 to n. Values are one less than the sum so
    that they can be used to index the elements of a vector storing
    the lower diagonal entries of a symmetric matrix that correspond
    to the diagonal entries
    """
    out = [int(0)] * n
    for i in range(n):
        out[i] = int(((i + 1) * (i + 2) / 2) - 1)
    return(out)
        
def ReadGRMBin(prefix, AllN = False):
    """
    read a GCTA binary GRM file storing relatedness values between individuals
    adapted from an R function on the GCTA website
    Davis McCarthy, February 2015
    """
    BinFileName  = prefix + ".grm.bin"
    NFileName = prefix + ".grm.N.bin"
    IDFileName = prefix + ".grm.id"
    dt = np.dtype('f4') # Relatedness is stored as a float of size 4 in the binary file
    entry_format = 'f' # N is stored as a float in the binary file
    entry_size = calcsize(entry_format)
    ## Read IDs
    ids = pd.read_csv(IDFileName, sep = '\t', header = None)
    ids_vec = ids.iloc[:,1]
    n = len(ids.index)
    ids_diag = ['NA' for x in range(n)]
    n_off = int((n * (n + 1) / 2) - n)
    ids_off = ['NA' for x in range(n_off)]
    ## Generate ids for relatedness values by concatenating individual IDs
    ticker = 0
    for i in range(n):
        for j in range(i):
            if i == j:
                ids_diag[i] = str(ids_vec[i])
            else:
                ids_off[ticker] = str(ids_vec[i]) + '_' + str(ids_vec[j])
                ticker += 1
    ## Read relatedness values
    grm = np.fromfile(BinFileName, dtype = dt)
    ## Read number of markers values
    if AllN:
        N = np.fromfile(NFileName, dtype = dt)
    else:
        with open(NFileName, mode='rb') as f:
            record = f.read(entry_size)
            N = unpack(entry_format, record)[0]
            N = int(N)
    i = sum_n_vec(n)
    out = {'diag': grm[i], 'off': np.delete(grm, i), 'id': ids, 'id_off': ids_off, 'id_diag': ids_diag, 'N': N}
    return(out)

def flatten(list_of_str):
    """
    flatten a list of list of strings into a list of the strings
    from http://stackoverflow.com/questions/5286541/how-can-i-flatten-lists-without-splitting-strings
    Davis McCarthy, June 2014
    """
    for x in list_of_str:
        if hasattr(x, '__iter__') and not isinstance(x, str):
            for y in flatten(x):
                yield y
        else:
            yield x


def plot_correlogram(df,figsize=(20,20)):
    """ 
    Create an n x n matrix of scatter plots for every
    combination of numeric columns in a pandas dataframe

    Credit: Corey Chivers (cjbayesian)
    https://github.com/cjbayesian
    Taken from:
    https://gist.github.com/cjbayesian/f0f127cc57f26d968f10
    2015-03-31
    """
    cols = list(df.columns[df.dtypes=='float64'])
    n = len(cols)
    fig, ax = plt.subplots(n,n,figsize=figsize)
    for i,y in enumerate(cols):
        for j,x in enumerate(cols):
            if i != n-1:
                ax[i,j].xaxis.set_ticklabels([])
            if j != 0:
                    ax[i,j].yaxis.set_ticklabels([])
            if i != j:
                try:
                    tmp = df[[x,y]].copy()
                    tmp.dropna(inplace=True)
                    ax[i,j].plot(tmp[x].values,tmp[y].values,'.',markersize=0.5,alpha=0.5,color='black')
                except:
                    pass
            else:
                midx = df[x].min() + (df[x].max() - df[x].min())/2.0
                midy = df[y].min() + (df[y].max() - df[y].min())/2.0
                ax[i,j].text(midx, midy, y.replace(' ','\n'),
                             horizontalalignment='center',
                             verticalalignment='center')
                ax[i,j].set_ylim((df[y].min(),df[y].max()))
                ax[i,j].set_xlim((df[x].min(),df[x].max()))
                            
            
