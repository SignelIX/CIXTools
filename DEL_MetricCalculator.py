import pandas as pd
import math
from scipy.stats import stats
from decimal import Decimal


########SELECTION NOT USED YET IN COUNTS##############
########only calculates trisynthons at the moment########
def Get_DELPopulationsStats ():
    def round_sig(x, sig=2):
        return round(x, sig - int(math.floor(math.log10(abs(x)))) - 1)
    tagdata = pd.read_csv('/Users/eric/Consulting/Clients/Cephalogix/UCB/DEL011 BBs and tags.csv')
    cyc_tagcts = tagdata.groupby(['Cycle']).agg({'Name': 'count'})
    BB_encodings = tagdata.groupby(['Cycle', 'Smiles']).agg({'Name': 'count'})

    seldata = pd.read_csv('/Users/eric/Consulting/Clients/Cephalogix/UCB/del_11_19_smiles.csv')
    libselcts = seldata.groupby(['library']).agg({'counts': 'sum'})
    libselct = libselcts.loc ['DEL11']['counts']

    n = libselct
    cyc1_N = cyc_tagcts.loc ['cyc1']['Name']
    cyc2_N = cyc_tagcts.loc ['cyc2']['Name']
    cyc3_N = cyc_tagcts.loc ['cyc3']['Name']
    N = cyc1_N * cyc2_N * cyc3_N
    recs = seldata.iloc [0]
    print (BB_encodings.loc ['cyc1'])
    for idx, recs in seldata.iterrows():
        library = recs['library']
        bb1_smiles = recs ['bb1_smiles']
        bb2_smiles = recs ['bb2_smiles']
        bb3_smiles = recs ['bb3_smiles']
        k = recs ['counts']
        bb1tagct = BB_encodings.loc ['cyc1'].loc[bb1_smiles]['Name']
        bb2tagct = BB_encodings.loc['cyc2'].loc[bb2_smiles]['Name']
        bb3tagct = BB_encodings.loc['cyc3'].loc[bb3_smiles]['Name']
        tagct = bb1tagct*bb2tagct*bb3tagct
        Kuai_LowCt = math.pow(math.sqrt (k+1)-2,2) * (N/tagct)/n
        Kuai_HighCt = math.pow(math.sqrt (k+1)+2,2) * (N/tagct)/n
        #Faver z score metric
        p_o = k/n
        p_i= tagct/N
        Faver_z_n = math.sqrt (p_i/(1-p_i))*((p_o/p_i)-1)
       # Faver confidence
        #z_alpha = 1.96  #95% CI, can adjust for different confidence level
        z_alpha =2.576 #98% CI, can adjust for different confidence level
        n_prime = n + math.pow (z_alpha,2)
        conf_p_o = (1 / n_prime) * (k + pow(z_alpha / 2, 2))
        Faver_low_conf = conf_p_o - z_alpha * math.sqrt((conf_p_o / n_prime) * (1 - conf_p_o))
        Faver_high_conf = conf_p_o + z_alpha * math.sqrt((conf_p_o / n_prime) * (1 - conf_p_o))
        Faver_z_n_low = math.sqrt (p_i/(1-p_i))*((Faver_low_conf/p_i)-1)
        Faver_z_n_high = math.sqrt(p_i / (1 - p_i)) * ((Faver_high_conf / p_i) - 1)
        Faver_enr_low = Faver_low_conf/p_i
        Faver_enr = p_o/p_i
        Faver_enr2 = conf_p_o / p_i
        Faver_enr_high = Faver_high_conf/p_i

        print ('General')

        print ('counts:', k)
        print('additional data [sel size, enc tags, population]', n, tagct, N)
        print ('sel norm counts: ', round_sig(k/n, 2))
        E = (tagct*(n/N))
        print('count_ratio: ', round((k-E)/E,2))
        print ('Enr, log Enr:', Faver_enr,  round (math.log10(Faver_enr), 2))
        crit_val =  stats.distributions.binom.ppf (1-(.05/n), n, tagct/N)
        print('CBV_ratio:')
        print (k/crit_val)
        print ('Kuai')
        print ( 'CI_enr: [', math.trunc (Kuai_LowCt), math.trunc(Kuai_HighCt) , ']')
        print ('Faver')
        print (  'CI_zscore: [', round(Faver_z_n_low,2),round(Faver_z_n,2), round(Faver_z_n_high,2), ']',
                 'CI_enr: [', round(Faver_enr_low,2), round(Faver_enr,2), round(Faver_enr_high,2), ']',
                 'CI_p: [', round_sig(Faver_low_conf,2),round_sig(p_o,2),  round_sig(Faver_high_conf,2), ']')

        h = 2*(math.asin(p_o)-math.asin(p_i))
        print ('Cohen\'s h:', h)
    #Kuai Metric
    # Kuai_EnrL = ((sqrt (k+1)-1)^2 * N/#encodings)/n
    # Kuai_EnrU = ((sqrt(k + 1) + 1) ^ 2 * N) / n
    # k = counts
    # N = # of total encodings (# of specific encodings)
    # n = selection size

    #Faver Metrics:
    #z_n = sqrt(p_i/(1-p_i))*((p_0/p_i)-1
    #p_i = expected freq, #of encodings/population of encodings =  tagct/N
    #p_0 = counts for BB combination/sample counts = k/n
    #Faver, Agresti-Coull uncertainty
    #p_0 +/- z_alpha * sqrt ((p_0/n_prime)(1-p_0)
    #p_0 = (1/n_prime)*(c_0 + pow(z_alpha/2,2)) c_0 equals k, obs counts
    #n_prime = n + pow(z_alpha,2)

    return




# Quantitative Comparison of Enrichment from DNA-Encoded Chemical Library Selections
# John C. Faver, Kevin Riehle, David R. Lancia, Jared B. J. Milbank, Christopher S. Kollmann, Nicholas Simmons, Zhifeng Yu, and Martin M. Matzuk
# ACS Combinatorial Science 2019 21 (2), 75-82
# DOI: 10.1021/acscombsci.8b00116
#Kuai L, O’Keeffe T, Arico-Muendel C. Randomness in DNA Encoded Library Selection Data Can Be Modeled for More Reliable Enrichment Calculation. SLAS DISCOVERY: Advancing the Science of Drug Discovery. 2018;23(5):405-416. doi:10.1177/2472555218757718