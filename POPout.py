#!/usr/bin/python3

import matplotlib 
import matplotlib.pyplot as plt
import sys
from collections import defaultdict as dd
import numpy as np
import scipy.stats as stats
import statsmodels.api as sm
import random
import math
import numpy as np
import pandas as pd
from scipy.stats import norm
from pathlib import Path
from collections import defaultdict as dd 

#import matplotlib
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes      
matplotlib.rcParams['xtick.labelsize'] = 24                                                                                                                                                                                 
matplotlib.rcParams['ytick.labelsize'] = 24                                                                                                                                                                                 
#matplotlib.rcParams['xtick.major.size'] = 18                                                                                                                                                                               
#matplotlib.rcParams['ytick.major.size'] = 18                                                                                                                                                                               
matplotlib.rcParams['axes.linewidth'] = 1.85  


plt.rc('text', usetex=True)
#matplotlib.use("Cairo")
import warnings
warnings.filterwarnings("ignore")


#############################################################################################
##################################      FUNCTIONS         ###################################
#############################################################################################


# ----------------------------------------
# Munge Functions (translations of R fns)
# ----------------------------------------

def est_prop_in_tail(effect_size, beta, r2, tail=0.01):
    """
    effect_size: length-2 array-like [lower_effect, upper_effect]
    beta: assumed effect size of rare large-effect variants
    r2: variance explained (from comboResult, first number)
    tail: tail proportion (default 1%)
    Returns: np.array([prop_lower, prop_upper])
    """
    kappa = norm.ppf(1 - tail)
    percentile_pushed_to_tail = norm.cdf(beta - kappa)

    #print(kappa) 
    #print(percentile_pushed_to_tail, 'hmmm')  
    # qnorm(percentile.pushed.to.tail, lower.tail=FALSE) -> norm.ppf(1 - p)
    Z1 = norm.ppf(1 - percentile_pushed_to_tail)
    Z0 = norm.ppf(1 - tail)

    # m.Y99.rare <- sqrt(r2) * dnorm(Z1) / percentile.pushed.to.tail
    mY99_rare = math.sqrt(r2) * norm.pdf(Z1) / max(percentile_pushed_to_tail, 1e-300)
    # m.Y99 <- sqrt(r2) * dnorm(Z0) / tail
    mY99 = math.sqrt(r2) * norm.pdf(Z0) / tail

    out = [0.0, 0.0]
    for i in range(2):
        if effect_size[i] > 0:
            denom = (mY99 - mY99_rare)
            
            #print(effect_size[i], mY99, mY99_rare) 
            #print(mY99 - mY99_rare) 

            if denom <= 0:
                out[i] = 1.0  # avoid negatives; conservative cap
            else:
                out[i] = effect_size[i] / denom

    return np.array(out, dtype=float)


def h2_rare_big(prop_in_tail, beta, tail=0.01, rare_maf=1e-4):
    """
    prop_in_tail: length-2 array-like, from est_prop_in_tail()
    beta: presumed effect size of rare large-effect variants
    tail: tail proportion
    rare_maf: rare MAF (R default 1e-4; the main loop calls with 1e-5)
    Returns dict matching R names:
      'h2','sum.rare.freq1','sum.rare.freq2','m1','m2','percentile.pushed.to.tail'
    """
    kappa = norm.ppf(1 - tail)
    percentile_pushed_to_tail = norm.cdf(beta - kappa)

    sum_rare_freq = []
    m = []
    var_rare = []
    for i in range(2):
        p = prop_in_tail[i]
        num = p * tail
        denom = p * tail + (1 - p) * percentile_pushed_to_tail
        denom = max(denom, 1e-300)
        srf = num / denom
        sum_rare_freq.append(srf)
        m_i = 0.5 * srf / rare_maf
        m.append(m_i)
        var_rare.append(m_i * 2 * rare_maf * (1 - rare_maf) * (beta ** 2))

    h2 = sum(var_rare) / (1 + sum(var_rare))
    return {
        'h2': h2,
        'sum.rare.freq1': sum_rare_freq[0],
        'sum.rare.freq2': sum_rare_freq[1],
        'm1': m[0],
        'm2': m[1],
        'percentile.pushed.to.tail': percentile_pushed_to_tail
    }


# Utility to parse a comma-sep numeric string and take the Nth (1-based like R)
def parse_comma_num(s, idx_1based):
    if pd.isna(s):
        return np.nan
    parts = str(s).split(',')
    if len(parts) < idx_1based:
        return np.nan
    try:
        return float(parts[idx_1based - 1])
    except Exception:
        return np.nan


# Summaries like in R for each beta
def summarize_beta(bidx, traits, x, rare_h2, rare_lower, rare_upper, theta_lower, theta_upper, p_tail_lower, p_tail_upper): 
    rows = []
    for i, trait in enumerate(traits):
        rows_trait = x[x.iloc[:, 0].astype(str) == trait]
        if len(rows_trait) != 9:
            rows.append([np.nan] * 13)
            continue

        h2_vals = rare_h2[i, bidx, :]
        th_lo = theta_lower[i, bidx, :]
        th_up = theta_upper[i, bidx, :]
        m_lo = rare_lower[i, bidx, :]
        m_up = rare_upper[i, bidx, :]
        p_lo_raw = p_tail_lower[i, bidx, :]
        p_up_raw = p_tail_upper[i, bidx, :]

        def qtiles(a, SHOW=False):
            
            a = a[~np.isnan(a)]
            if a.size == 0:
                return [np.nan, np.nan, np.nan]
            
            if SHOW: 
                print('clive', len(a), np.mean(a)) 
            return np.round(np.quantile(a, [0.5, 0.025, 0.975]), 3).tolist()

        def med(a):
            a = a[~np.isnan(a)]
            return float(np.round(np.median(a), 1)) if a.size else np.nan

        def mean_bool(cond):
            cond = cond[~np.isnan(cond)]
            return float(np.round(np.mean(cond.astype(float)), 2)) if cond.size else np.nan

        row = []
        row += qtiles(h2_vals, True)
        row += qtiles(th_lo)
        row += qtiles(th_up)
        row += [med(m_lo), med(m_up)]
        row += [mean_bool(p_lo_raw < 1), mean_bool(p_up_raw < 1)]
        

        rows.append(row)

    df = pd.DataFrame(rows, columns=[
        'h2','h2.2.5%','h2.97.5%',
        'theta.lt','theta.lt.2.5%','theta.lt.97.5%',
        'theta.ut','theta.ut.2.5%','theta.ut.97.5%',
        'm.lt','m.ut','val.lt','val.ut'
    ])
    df.insert(0, 'trait', traits)
    return df


def get_quantiles(X): 
    qx = [x for x in X if not np.isnan(x)] 
    if len(qx) == 0: return [np.nan, np.nan, np.nan]
    return np.round(np.quantile(qx, [0.5, 0.025, 0.975]), 3).tolist()

def cond_bool(X): 
    qx = [x for x in X if not np.isnan(x)]
    return round(np.mean([float(x < 1) for x in qx]),3) 




#############################################################################################
##################################       CLASSES          ###################################
#############################################################################################




#############################################################################################
##################################  PROGRESS STREAM  ########################################
#############################################################################################


class PopProgress:
    def __init__(self, args, command_line, MODE='Analyze Tails'): 
        self.args = args
        if args.silent: self.ACTIVE = False 
        else:           self.ACTIVE = True 
        self.out2  = sys.stderr 
        self.space = '' 
        self.show('\nPopOut Begins:  '+command_line+'\n')
        self.show('         Mode:  '+MODE+'\n') 
        if len(args.popfiles) > 0: self.show('  Input Files:  '+",".join([sf.name.split('/')[-1] for sf in args.popfiles])+'\n') 
        else: self.show('  Input Files:  '+",".join([sf.name.split('/')[-1] for sf in args.makePOPfile])+'\n') 
        self.show('Output Prefix:  '+self.args.out+'\n\n')
        self.loc = None
        self.spl = '' 

    def initialize(self, f_name, t, i): 
        if self.loc is None: 
            self.show('Initializing Trait Analysis: '+t+'\n') 
        else: 
            self.show('Complete\nInitializing Trait Analysis: '+t+'\n') 

        self.space = '         '
        self.show('Reading Input File: '+f_name+'\n') 
        


    def begin_analysis(self, tail = 1, COMMAND='NA'): 
        if tail >= 1: ts = str(int(tail))+'%' 
        else:         ts = str(tail)+'%'  
        if self.loc is not None: self.show('Complete\n') 
        self.space = '         ' 
        if COMMAND == 'PLOT': self.show('Collating Trait Result & Drawing Plot......') 
        elif COMMAND == 'PREDS': self.show('Evaluating PRS Tail Performance.......') 
        elif COMMAND == 'MUNGE': self.show('Estimating High Effect Rare Heritability...') 
        else:    self.show('Beginning '+ts+' Tail Analysis.....') 
        self.loc, self.space = COMMAND, '' 

    def show(self, msg, space=''):
        if space == 'NA': myspace = '' 
        else:             myspace = self.space
        if self.ACTIVE: 
            self.out2.write(myspace+msg)
            self.out2.flush() 

    def step(self, dots): self.show(''.join(['.' for x in range(dots)])) 
    
    def complete_step(self): self.show('...Complete\n') 

    def complete_analysis(self): 
        if self.loc is not None: self.show('...Finished\n','NA') 
        self.show('\nPOPout Completed\n','NA') 



#############################################################################################
##################################     TRAIT TABLE        ###################################
#############################################################################################


class TailTable:
    def __init__(self, ax,  result, clr='white', PETT=False):
        self.ax, self.result, self.rows, self.cols, self.clr  = ax, result, [], [], clr

    def get_loc(self,X,Y):
        x1,x2 = X[0]/100.0 , X[1]/100.0
        y1,y2 = Y[0]/100.0 , Y[1]/100.0
        return (x1,y1,x2-x1,y2-y1)

    def add_row(self,row_data,X=None,Y=None,COLORS=[],WIDTHS=[],FS=12,ALPHA=0,TITLE=False, CL = 'center',CLEAR=False):
        if X == None: X = self.X_SPAN
        if Y == None: Y = self.Y_SPAN
        cL,rL,rD,xL = CL,None,[row_data],len(row_data)
        bl = self.get_loc(X,Y)
        while len(WIDTHS) < len(row_data): WIDTHS.append(10)
        while len(COLORS) < len(row_data): COLORS.append('white')
        if CL != 'center': row = self.ax.table(cellText=rD,rowLabels=rL,cellColours = [COLORS[0:xL]],colWidths=WIDTHS[0:xL],bbox=bl,loc = cL, cellLoc=cL, alpha = ALPHA, clip_on=False) 
        else:              row = self.ax.table(cellText=rD,rowLabels=rL,cellColours = [COLORS[0:xL]],colWidths=WIDTHS[0:xL],bbox=bl,cellLoc=cL, alpha=ALPHA, clip_on=False) 
        row.auto_set_font_size(False)
        row.set_fontsize(FS)
        table_props = row.properties()
        if ALPHA > 0: 
            for cell in row._cells: row._cells[cell].set_alpha(ALPHA)
        self.rows.append(row)



    def make_popout(self,SK, c1 = 'lightgrey', c2='snow'):
        self.widths = [16,12,35,19,18]
        self.vwid = [0] + [sum(self.widths[0:i+1]) for i in range(len(self.widths))]
        self.y2, fs1, fs2, fs3 = 70, 30, 20, 18
        header = ['Tail','Size','POPout Effect (CI)','P-Value','Odds Ratio']
        for i,h in enumerate(header):
            if i in [2,4]: self.add_row([h], COLORS=[c1], X=(self.vwid[i],self.vwid[i+1]), Y=(self.y2,100), FS=fs3, WIDTHS=[self.widths[i]], TITLE=True)
            else:    self.add_row([h], COLORS=[c1], X=(self.vwid[i],self.vwid[i+1]), Y=(self.y2,100), FS=fs2, WIDTHS=[self.widths[i]], TITLE=True)
        sk = SK*100 
        if sk == int(sk): sk_str = str(int(sk))+'\%' 
        else:             sk_str = str(sk)+'\%' 
        for i,(name,loc) in enumerate([['upper',1],['lower',0]]): 
            if loc == 1: ef,e1,e2 = [str(round(-1*x,2)) for x in self.result['effects'][loc]]
            else:        ef,e1,e2 = [str(round(x,2)) for x in self.result['effects'][loc]]
            
            if e1 < e2: e_str = ef + ' ('+e1+','+e2+')'
            else:       e_str = ef + ' ('+e2+','+e1+')'

            p1 = self.result['pvals'][loc]
            if p1 < 0.001: p1 = str('%5.2e' % p1) 
            else:          p1 = str(round(p1,3)) 
            if i == 0: p2 = round(self.result['preds']['b2'].upperTailOdds,2) 
            else:      p2 = round(self.result['preds']['b1'].lowerTailOdds,2) 
            my_data = [name, sk_str,e_str,p1,str(p2)] 
            self.add_row(my_data, COLORS=[c2 for c in my_data], X=(self.vwid[0],self.vwid[5]), Y=(self.y2-35,self.y2), FS=fs3, WIDTHS=self.widths[0:5], TITLE=True)
            self.y2 -= 35
        self.ax.axis('off') 


    def make_summary(self,SK, c1 = 'lightgrey', c2='snow', fs1=30, fs2=20, fs3 = 18):
        R = self.result
        self.widths, self.vwid, self.y2 = [50,50], [0,100], 100 
        step = 17 
        self.add_row(['Trait Summary'], COLORS=[c1], X=(self.vwid[0],self.vwid[1]), Y=(self.y2-step,self.y2), FS=fs3, WIDTHS=[100], TITLE=True)
        
        r95,p95 = R['R'][95] 
        r5,p5   = R['R'][5]
        
        #Collating Trait Result & Drawing Plot......dict_keys(['tot', 25, 10, 5, 1, 50, 75, 90, 95, 99])


        #for i,(name,v) in enumerate([['Samples',R['len']],['$R^2$',R['r2']],['QC',R['QC']],['$h^2_{HER}$','h2her']]): 
        for i,(name,v) in enumerate([['Samples',R['len']],['$R^2$',R['r2']],['QC',R['QC']],['$R_{<5\%}$',R['R'][5]],['$R_{>95\%}$',R['R'][95]]]): 
            self.y2 -= step
            if i == 1: my_row = [name, round(v,2)]
            elif i < 3: my_row = [name, v] 
            
            else: 
                my_row = [name, 'NA'] 
                if v == 'h2her': 
                    try: 
                        beta = sorted([k for k in R['munge'].estimates.keys()])[-1]  
                        val = R['munge'].estimates[beta]['h2'][0] 
                        my_row = [name+'\n$(\\beta='+str(beta)+')$', val] 
                    except: pass
                else: 
                    r = str(round(v[0],2)) 
                    if v[1] > 0.05:  r = 'NS' 
                    my_row = [name, r] 
            self.add_row(my_row, COLORS=[c2,c2], X=(self.vwid[0],self.vwid[1]), Y=(self.y2-step,self.y2), FS=fs3, WIDTHS=self.widths, TITLE=True)
        self.ax.axis('off')

        return





#############################################################################################
##################################     TRAIT PLOT         ###################################
#############################################################################################



class TailPlot:
    def __init__(self, args, progress): 
        self.args, self.progress =  args, progress
        self.fs1, self.fs2, self.fs3 = 34, 24, 15

        #print(len(self.args.popfiles)) 


    def setup(self):
        self.WD, self.HT, self.rows, self.cols = 10, 12, 3, 2  
        self.fig = matplotlib.pyplot.gcf()
        self.fig.set_size_inches(self.WD, self.HT) 
        self.axes, self.ax_index = [], 0  
        self.rows, self.cols = 16, 7 
        self.axes.append(plt.subplot2grid((self.rows, self.cols), (0,0),     rowspan=5, colspan=4)) 
        self.axes.append(plt.subplot2grid((self.rows, self.cols), (0,5),     rowspan=5, colspan=2)) 
        self.axes.append(plt.subplot2grid((self.rows, self.cols), (7,0),     rowspan=2, colspan=7)) 
        self.axes.append(plt.subplot2grid((self.rows, self.cols), (11,0),     rowspan=3, colspan=3)) 
        self.axes.append(plt.subplot2grid((self.rows, self.cols), (11,4),     rowspan=3, colspan=3)) 
        self.axes.append(plt.subplot2grid((self.rows, self.cols), (14,0),     rowspan=2, colspan=3)) 
        self.axes.append(plt.subplot2grid((self.rows, self.cols), (14,4),     rowspan=2, colspan=3)) 
        return

    def create_trait_plot(self, pop, name): 
        
        self.pop, self.name = pop, name
        self.setup() 
        self.draw_curve(self.axes[0], pop.results[name], name) 
        self.draw_tables(self.axes[1], self.axes[2], pop.results[name], name)   
        self.draw_preds(self.axes[3::], pop.results[name], name) 
        self.finish(name) 
        return self


    def draw_tables(self, ax1, ax2, result, name): 
        popTable = TailTable(ax1, result).make_summary(self.pop.SK) 
        popTable = TailTable(ax2, result).make_popout(self.pop.SK) 
        ax1.set_xticks([]) 
        ax2.set_xticks([]) 


    def get_expected_pred(self, R2, t): 
        if R2==1:
            if t==25: return [1.4,1.3,1.2,1.1,1.1,1.0,1.0,0.9,0.9,0.8,0.7],[0.7,0.8,0.9,0.9,1.0,1.0,1.0,1.1,1.2,1.2,1.4]                                                                              
            elif t==10: return [1.5,1.3,1.2,1.1,1.0,1.0,1.0,0.9,0.8,0.8,0.7],[0.7,0.8,0.8,0.9,1.0,1.0,1.1,1.1,1.2,1.3,1.5]
            elif t==5: return [1.5,1.3,1.2,1.1,1.1,1.0,1.0,0.9,0.8,0.7,0.6],[0.7,0.7,0.8,0.9,1.0,1.0,1.1,1.1,1.2,1.3,1.6]
            elif t==1: return [1.8,1.5,1.3,1.1,1.1,1.0,1.0,0.9,0.9,0.7,0.6],[0.6,0.6,0.7,0.8,1.0,1.0,1.1,1.1,1.2,1.4,1.7]
            elif t==0.5: return [1.8,1.5,1.3,1.1,1.1,1.0,1.0,1.0,0.8,0.7,0.5],[0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.1,1.2,1.4,1.7]
        elif R2==2:
           if t==25: return [1.6,1.4,1.2,1.1,1.1,1.0,0.9,0.9,0.8,0.7,0.6],[0.6,0.7,0.8,0.9,0.9,1.0,1.1,1.1,1.2,1.4,1.6]
           elif t==10: return [1.7,1.4,1.3,1.2,1.1,1.0,0.9,0.9,0.8,0.7,0.6],[0.5,0.7,0.8,0.9,0.9,1.0,1.1,1.2,1.3,1.4,1.7]
           elif t==5: return [1.9,1.5,1.3,1.2,1.1,1.0,1.0,0.9,0.8,0.7,0.5],[0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.5,1.9]
           elif t==1: return [2.08,1.6,1.3,1.3,1.1,1.0,0.9,0.8,0.7,0.6,0.5],[0.4,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,2.13]
           elif t==0.5: return [2.23,1.7,1.4,1.3,1.2,1.0,0.9,0.8,0.7,0.6,0.5],[0.4,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.7,2.21]
        elif R2==3:
           if t==25: return [1.8,1.5,1.3,1.2,1.1,1.0,0.9,0.8,0.8,0.7,0.5],[0.5,0.7,0.8,0.9,0.9,1.0,1.1,1.2,1.3,1.5,1.8]
           elif t==10: return [2.0,1.5,1.3,1.2,1.1,1.0,0.9,0.8,0.7,0.6,0.5],[0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.5,2.0]
           elif t==5: return [2.13,1.6,1.3,1.2,1.1,1.0,0.9,0.8,0.7,0.6,0.4],[0.4,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,2.13]
           elif t==1: return [2.53,1.8,1.5,1.3,1.2,1.0,0.9,0.8,0.6,0.5,0.4],[0.4,0.5,0.7,0.8,0.9,1.0,1.1,1.3,1.5,1.9,2.65]
           elif t==0.5: return [2.78,2.0,1.6,1.3,1.1,1.0,0.9,0.8,0.7,0.5,0.3],[0.4,0.5,0.7,0.8,0.9,1.0,1.1,1.4,1.5,2.0,2.81]
        elif R2==4:
           if t==25: return [2.02,1.5,1.3,1.2,1.1,1.0,0.9,0.8,0.7,0.6,0.5],[0.5,0.6,0.8,0.8,0.9,1.0,1.1,1.2,1.3,1.6,2.02]
           elif t==10: return [2.16,1.6,1.4,1.2,1.1,1.0,0.9,0.8,0.7,0.6,0.4],[0.4,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.7,2.17]
           elif t==5: return [2.4,1.7,1.4,1.2,1.1,1.0,0.8,0.8,0.7,0.5,0.4],[0.4,0.5,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.7,2.34]
           elif t==1: return [3.04,2.05,1.5,1.3,1.2,1.0,0.8,0.7,0.6,0.5,0.3],[0.3,0.5,0.6,0.8,1.0,1.0,1.2,1.3,1.6,2.07,3.0]
           elif t==0.5: return [3.47,2.22,1.5,1.3,1.2,1.0,0.9,0.7,0.6,0.5,0.3],[0.3,0.5,0.6,0.9,1.0,1.0,1.2,1.3,1.8,2.14,3.34]
        elif R2==5:
           if t==25: return [2.22,1.7,1.4,1.2,1.1,1.0,0.9,0.8,0.7,0.6,0.4],[0.4,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,2.16]
           elif t==10: return [2.45,1.8,1.5,1.3,1.1,1.0,0.9,0.8,0.7,0.5,0.4],[0.4,0.5,0.7,0.8,0.9,1.0,1.1,1.3,1.5,1.7,2.4]
           elif t==5: return [2.74,1.9,1.5,1.3,1.1,1.0,0.9,0.8,0.6,0.5,0.3],[0.3,0.5,0.6,0.8,0.9,1.0,1.1,1.3,1.5,1.9,2.65]
           elif t==1: return [3.59,2.25,1.7,1.4,1.2,1.0,0.9,0.8,0.6,0.4,0.3],[0.3,0.4,0.6,0.7,0.8,1.0,1.2,1.3,1.6,2.13,3.25]
           elif t==0.5: return [4.01,2.34,1.8,1.4,1.2,1.0,0.8,0.7,0.5,0.4,0.2],[0.2,0.4,0.5,0.7,0.8,1.0,1.2,1.4,1.7,2.24,3.47]
        elif R2==6:
           if t==25: return [2.36,1.7,1.4,1.3,1.1,1.0,0.9,0.8,0.7,0.6,0.4],[0.4,0.6,0.7,0.8,0.9,1.0,1.1,1.3,1.4,1.7,2.38]
           elif t==10: return [2.67,1.9,1.5,1.3,1.1,1.0,0.9,0.8,0.7,0.5,0.3],[0.4,0.5,0.7,0.8,0.9,1.0,1.1,1.3,1.5,1.9,2.68]
           elif t==5: return [3.0,2.07,1.6,1.4,1.2,1.0,0.9,0.8,0.6,0.5,0.3],[0.3,0.5,0.6,0.7,0.8,1.0,1.2,1.4,1.6,2.03,3.02]
           elif t==1: return [3.86,2.4,1.8,1.6,1.2,1.0,0.8,0.7,0.5,0.4,0.2],[0.2,0.4,0.5,0.7,0.8,1.0,1.3,1.5,1.7,2.41,4.07]
           elif t==0.5: return [4.08,2.45,1.8,1.6,1.1,1.0,0.8,0.6,0.5,0.4,0.2],[0.2,0.4,0.6,0.6,0.8,1.0,1.3,1.5,1.8,2.51,4.65]
        elif R2==7:
           if t==25: return [2.53,1.8,1.5,1.3,1.1,1.0,0.9,0.8,0.7,0.5,0.4],[0.4,0.5,0.7,0.8,0.9,1.0,1.1,1.3,1.5,1.8,2.52]
           elif t==10: return [2.89,2.0,1.6,1.3,1.1,1.0,0.9,0.7,0.6,0.5,0.3],[0.3,0.5,0.6,0.8,0.9,1.0,1.1,1.3,1.6,1.9,2.86]
           elif t==5: return [3.25,2.12,1.6,1.4,1.2,1.0,0.9,0.7,0.6,0.4,0.3],[0.3,0.4,0.6,0.7,0.9,1.0,1.2,1.4,1.6,2.1,3.21]
           elif t==1: return [4.52,2.73,2.0,1.5,1.3,1.0,0.9,0.7,0.5,0.4,0.2],[0.2,0.3,0.5,0.7,0.8,1.0,1.2,1.4,1.9,2.46,4.29]
           elif t==0.5: return [5.05,2.9,2.0,1.6,1.3,1.0,0.8,0.7,0.5,0.4,0.2],[0.2,0.3,0.5,0.7,0.8,1.0,1.2,1.5,2.01,2.71,4.97]
        elif R2==8:
           if t==25: return [2.69,1.9,1.5,1.3,1.1,1.0,0.9,0.8,0.6,0.5,0.3],[0.3,0.5,0.6,0.8,0.9,1.0,1.1,1.3,1.5,1.9,2.71]
           elif t==10: return [3.03,2.03,1.6,1.3,1.1,1.0,0.9,0.7,0.6,0.5,0.3],[0.3,0.4,0.6,0.7,0.9,1.0,1.1,1.3,1.6,2.06,3.04]
           elif t==5: return [3.41,2.2,1.7,1.4,1.1,1.0,0.8,0.7,0.6,0.4,0.2],[0.2,0.4,0.5,0.7,0.8,1.0,1.1,1.4,1.7,2.2,3.38]
           elif t==1: return [4.74,2.88,2.0,1.6,1.2,1.0,0.8,0.7,0.6,0.4,0.2],[0.1,0.3,0.5,0.7,0.8,1.0,1.2,1.6,2.08,2.82,4.86]
           elif t==0.5: return [5.35,3.2,2.0,1.7,1.3,1.0,0.8,0.7,0.5,0.3,0.2],[0.1,0.3,0.5,0.7,0.9,1.0,1.3,1.7,2.4,3.36,6.01]
        elif R2==9:
           if t==25: return [2.9,2.0,1.6,1.3,1.1,1.0,0.9,0.7,0.6,0.5,0.3],[0.3,0.5,0.6,0.7,0.9,1.0,1.2,1.3,1.6,2.0,2.91]
           elif t==10: return [3.26,2.15,1.7,1.4,1.2,1.0,0.8,0.7,0.6,0.4,0.3],[0.3,0.4,0.6,0.7,0.8,1.0,1.2,1.4,1.6,2.16,3.27]
           elif t==5: return [3.62,2.28,1.7,1.4,1.2,1.0,0.8,0.7,0.5,0.4,0.2],[0.2,0.4,0.5,0.7,0.8,1.0,1.2,1.4,1.7,2.32,3.7]
           elif t==1: return [4.81,2.69,1.9,1.5,1.2,1.0,0.7,0.6,0.4,0.3,0.2],[0.2,0.3,0.4,0.6,0.8,1.0,1.2,1.5,1.9,2.8,5.09]
           elif t==0.5: return [5.49,2.94,2.0,1.6,1.2,1.0,0.7,0.6,0.4,0.2,0.1],[0.1,0.3,0.5,0.6,0.8,1.0,1.2,1.6,2.17,3.11,6.1]
        elif R2==10:
           if t==25: return [3.09,2.04,1.6,1.3,1.2,1.0,0.9,0.7,0.6,0.5,0.3],[0.3,0.5,0.6,0.7,0.9,1.0,1.2,1.4,1.6,2.05,3.06]
           elif t==10: return [3.54,2.26,1.7,1.4,1.2,1.0,0.9,0.7,0.6,0.4,0.2],[0.2,0.4,0.6,0.7,0.8,1.0,1.2,1.4,1.7,2.27,3.54]
           elif t==5: return [4.05,2.47,1.8,1.4,1.2,1.0,0.8,0.6,0.5,0.4,0.2],[0.2,0.4,0.5,0.7,0.8,1.0,1.2,1.5,1.8,2.48,4.01]
           elif t==1: return [5.85,3.17,2.1,1.6,1.3,1.0,0.8,0.6,0.4,0.3,0.2],[0.1,0.3,0.4,0.6,0.8,1.0,1.2,1.6,2.11,3.05,5.7]
           elif t==0.5: return [7.2,3.85,2.47,1.8,1.3,1.0,0.9,0.6,0.4,0.3,0.1],[0.1,0.3,0.4,0.6,0.8,1.0,1.3,1.7,2.25,3.38,6.58]
        elif R2==11:
           if t==25: return [3.27,2.14,1.6,1.4,1.2,1.0,0.9,0.7,0.6,0.4,0.3],[0.3,0.4,0.6,0.7,0.9,1.0,1.2,1.4,1.6,2.13,3.24]
           elif t==10: return [3.76,2.34,1.7,1.4,1.2,1.0,0.8,0.7,0.5,0.4,0.2],[0.2,0.4,0.5,0.7,0.8,1.0,1.2,1.4,1.8,2.39,3.76]
           elif t==5: return [4.33,2.59,1.8,1.5,1.2,1.0,0.8,0.7,0.5,0.3,0.2],[0.2,0.3,0.5,0.7,0.8,1.0,1.2,1.5,1.9,2.65,4.38]
           elif t==1: return [6.29,3.4,2.24,1.7,1.3,1.0,0.8,0.6,0.4,0.2,0.1],[0.1,0.2,0.4,0.6,0.8,1.0,1.3,1.7,2.18,3.44,6.72]
           elif t==0.5: return [7.58,3.91,2.59,1.9,1.4,1.0,0.8,0.6,0.5,0.2,0.1],[0.1,0.2,0.4,0.6,0.8,1.0,1.3,1.8,2.4,4.05,8.28]
        elif R2==12:
           if t==25: return [3.43,2.18,1.7,1.4,1.2,1.0,0.8,0.7,0.6,0.4,0.2],[0.2,0.4,0.6,0.7,0.9,1.0,1.2,1.4,1.7,2.24,3.45]
           elif t==10: return [4.0,2.44,1.8,1.5,1.2,1.0,0.8,0.7,0.5,0.4,0.2],[0.2,0.4,0.5,0.7,0.8,1.0,1.2,1.5,1.8,2.46,3.99]
           elif t==5: return [4.7,2.77,2.0,1.6,1.2,1.0,0.8,0.7,0.5,0.3,0.1],[0.2,0.3,0.5,0.6,0.8,1.0,1.2,1.6,1.9,2.72,4.63]
           elif t==1: return [7.25,3.61,2.4,1.7,1.3,1.0,0.8,0.6,0.4,0.2,0.1],[0.1,0.2,0.4,0.6,0.8,1.0,1.3,1.7,2.29,3.56,6.83]
           elif t==0.5: return [8.59,4.04,2.56,1.8,1.2,1.0,0.7,0.6,0.4,0.2,0.1],[0.1,0.2,0.4,0.5,0.8,1.0,1.3,1.8,2.54,4.09,8.46]
        elif R2==13:
           if t==25: return [3.67,2.3,1.7,1.4,1.2,1.0,0.8,0.7,0.6,0.4,0.2],[0.2,0.4,0.6,0.7,0.8,1.0,1.2,1.4,1.7,2.28,3.64]
           elif t==10: return [4.2,2.54,1.8,1.5,1.2,1.0,0.8,0.7,0.5,0.3,0.2],[0.2,0.3,0.5,0.7,0.8,1.0,1.2,1.5,1.9,2.56,4.25]
           elif t==5: return [4.88,2.8,2.0,1.5,1.2,1.0,0.8,0.6,0.4,0.3,0.1],[0.1,0.3,0.5,0.6,0.8,1.0,1.2,1.5,2.0,2.81,4.92]
           elif t==1: return [7.26,3.55,2.24,1.6,1.3,1.0,0.7,0.5,0.3,0.2,0.1],[0.1,0.2,0.4,0.5,0.7,1.0,1.3,1.7,2.29,3.61,7.39]
           elif t==0.5: return [8.88,4.02,2.48,1.7,1.4,1.0,0.7,0.5,0.3,0.2,0.1],[0.1,0.2,0.3,0.5,0.7,1.0,1.3,1.8,2.52,4.09,9.28]
        elif R2==14:
           if t==25: return [3.9,2.38,1.8,1.4,1.2,1.0,0.8,0.7,0.6,0.4,0.2],[0.2,0.4,0.5,0.7,0.8,1.0,1.2,1.4,1.7,2.34,3.81]
           elif t==10: return [4.58,2.64,1.9,1.5,1.2,1.0,0.8,0.6,0.5,0.3,0.2],[0.2,0.3,0.5,0.7,0.8,1.0,1.2,1.5,1.9,2.65,4.52]
           elif t==5: return [5.42,2.95,2.03,1.6,1.3,1.0,0.8,0.6,0.4,0.3,0.1],[0.1,0.3,0.5,0.6,0.8,1.0,1.3,1.6,2.12,3.03,5.39]
           elif t==1: return [8.32,3.89,2.5,1.8,1.4,1.0,0.7,0.5,0.3,0.2,0.1],[0.1,0.2,0.3,0.5,0.7,1.0,1.4,1.8,2.49,4.05,8.21]
           elif t==0.5: return [10.3,4.47,2.75,1.9,1.4,1.0,0.7,0.5,0.3,0.2,0.1],[0.1,0.2,0.3,0.5,0.7,1.0,1.4,2.0,2.83,4.55,10.12]
        elif R2==15:
           if t==25: return [4.15,2.48,1.8,1.5,1.2,1.0,0.8,0.7,0.5,0.4,0.2],[0.2,0.4,0.5,0.7,0.8,1.0,1.2,1.4,1.8,2.46,4.09]
           elif t==10: return [5.04,2.83,2.01,1.6,1.3,1.0,0.8,0.6,0.5,0.3,0.2],[0.1,0.3,0.5,0.6,0.8,1.0,1.3,1.6,2.0,2.79,4.91]
           elif t==5: return [5.88,3.17,2.19,1.6,1.3,1.0,0.8,0.6,0.5,0.3,0.1],[0.1,0.3,0.4,0.6,0.7,1.0,1.3,1.6,2.15,3.18,5.82]
           elif t==1: return [9.37,4.53,2.67,1.9,1.3,1.0,0.8,0.6,0.4,0.2,0.1],[0.1,0.2,0.3,0.5,0.7,1.0,1.3,1.8,2.53,4.23,9.27]
           elif t==0.5: return [12.13,4.92,3.05,1.9,1.3,1.0,0.7,0.5,0.3,0.2,0.0],[0.1,0.1,0.3,0.6,0.7,1.0,1.4,2.0,2.72,4.93,11.76]
        elif R2==16:
           if t==25: return [4.2,2.55,1.8,1.4,1.2,1.0,0.8,0.7,0.5,0.4,0.2],[0.2,0.4,0.5,0.7,0.8,1.0,1.2,1.5,1.8,2.54,4.3]
           elif t==10: return [5.07,2.9,2.0,1.6,1.3,1.0,0.8,0.6,0.4,0.3,0.1],[0.1,0.3,0.5,0.6,0.8,1.0,1.2,1.6,2.04,2.87,5.07]
           elif t==5: return [6.13,3.32,2.21,1.6,1.3,1.0,0.8,0.6,0.4,0.3,0.1],[0.1,0.2,0.4,0.6,0.8,1.0,1.3,1.7,2.26,3.3,6.12]
           elif t==1: return [9.38,4.4,2.57,1.8,1.4,1.0,0.7,0.5,0.3,0.2,0.1],[0.0,0.2,0.3,0.5,0.7,1.0,1.4,1.9,2.75,4.41,9.98]
           elif t==0.5: return [12.18,5.25,3.0,2.0,1.5,1.0,0.7,0.5,0.2,0.2,0.0],[0.0,0.1,0.3,0.5,0.7,1.0,1.5,2.09,3.07,5.33,13.29]
        elif R2==16:
           if t==25: return [4.2,2.55,1.8,1.4,1.2,1.0,0.8,0.7,0.5,0.4,0.2],[0.2,0.4,0.5,0.7,0.8,1.0,1.2,1.5,1.8,2.54,4.3]
           elif t==10: return [5.07,2.9,2.0,1.6,1.3,1.0,0.8,0.6,0.4,0.3,0.1],[0.1,0.3,0.5,0.6,0.8,1.0,1.2,1.6,2.04,2.87,5.07]
           elif t==5: return [6.13,3.32,2.21,1.6,1.3,1.0,0.8,0.6,0.4,0.3,0.1],[0.1,0.2,0.4,0.6,0.8,1.0,1.3,1.7,2.26,3.3,6.12]
           elif t==1: return [9.38,4.4,2.57,1.8,1.4,1.0,0.7,0.5,0.3,0.2,0.1],[0.0,0.2,0.3,0.5,0.7,1.0,1.4,1.9,2.75,4.41,9.98]
           elif t==0.5: return [12.18,5.25,3.0,2.0,1.5,1.0,0.7,0.5,0.2,0.2,0.0],[0.0,0.1,0.3,0.5,0.7,1.0,1.5,2.09,3.07,5.33,13.29]
        elif R2==17:
           if t==50: return [4.71,2.62,1.8,1.5,1.2,1.0,0.8,0.7,0.5,0.4,0.2],[0.2,0.4,0.5,0.7,0.8,1.0,1.2,1.4,1.8,2.6,4.51]
           elif t==25: return [4.56,2.61,1.9,1.5,1.2,1.0,0.8,0.7,0.5,0.3,0.2],[0.2,0.3,0.5,0.6,0.8,1.0,1.2,1.5,1.9,2.58,4.43]
           elif t==10: return [5.51,3.02,2.13,1.6,1.3,1.0,0.8,0.6,0.5,0.3,0.1],[0.1,0.3,0.4,0.6,0.8,1.0,1.3,1.6,2.07,2.99,5.3]
           elif t==5: return [6.59,3.39,2.32,1.7,1.3,1.0,0.8,0.6,0.4,0.2,0.1],[0.1,0.2,0.4,0.6,0.8,1.0,1.3,1.7,2.25,3.4,6.43]
           elif t==1: return [10.78,4.46,2.75,1.8,1.3,1.0,0.8,0.4,0.3,0.2,0.1],[0.1,0.2,0.3,0.5,0.7,1.0,1.5,2.02,2.91,4.76,11.22]
           elif t==0.5: return [13.3,4.92,3.09,1.9,1.3,1.0,0.7,0.4,0.3,0.2,0.0],[0.0,0.1,0.3,0.4,0.7,1.0,1.6,1.9,3.18,5.15,13.04]
        elif R2==18:
           if t==25: return [4.8,2.72,1.9,1.5,1.2,1.0,0.8,0.7,0.5,0.3,0.2],[0.2,0.3,0.5,0.7,0.8,1.0,1.2,1.5,1.9,2.72,4.86]
           elif t==10: return [5.7,3.1,2.13,1.6,1.3,1.0,0.8,0.6,0.4,0.3,0.1],[0.1,0.3,0.4,0.6,0.8,1.0,1.3,1.6,2.17,3.14,5.83]
           elif t==5: return [7.0,3.57,2.35,1.7,1.3,1.0,0.8,0.6,0.4,0.2,0.1],[0.1,0.2,0.4,0.6,0.8,1.0,1.3,1.7,2.38,3.59,7.1]
           elif t==1: return [11.75,5.07,2.98,2.0,1.4,1.0,0.7,0.5,0.3,0.1,0.0],[0.0,0.1,0.3,0.5,0.7,1.0,1.3,2.0,2.91,4.84,11.78]
           elif t==0.5: return [15.57,6.08,3.43,2.29,1.6,1.0,0.8,0.5,0.2,0.1,0.0],[0.0,0.1,0.2,0.4,0.6,1.0,1.3,2.09,3.12,5.58,14.67]
        elif R2==19:
           if t==25: return [5.08,2.83,2.0,1.5,1.2,1.0,0.8,0.6,0.5,0.3,0.2],[0.1,0.3,0.5,0.6,0.8,1.0,1.2,1.5,2.0,2.79,5.03]
           elif t==10: return [6.07,3.25,2.2,1.7,1.3,1.0,0.8,0.6,0.4,0.3,0.1],[0.1,0.2,0.4,0.6,0.8,1.0,1.3,1.6,2.17,3.23,6.04]
           elif t==5: return [7.52,3.8,2.46,1.8,1.3,1.0,0.8,0.6,0.4,0.2,0.1],[0.1,0.2,0.4,0.5,0.7,1.0,1.3,1.8,2.46,3.83,7.47]
           elif t==1: return [13.24,5.59,3.11,2.14,1.5,1.0,0.7,0.4,0.3,0.1,0.0],[0.1,0.1,0.3,0.5,0.7,1.0,1.5,2.14,3.13,5.5,13.01]
           elif t==0.5: return [17.66,6.59,3.58,2.34,1.5,1.0,0.7,0.4,0.3,0.1,0.0],[0.1,0.1,0.3,0.4,0.7,1.0,1.7,2.43,3.63,6.61,17.63]
        elif R2==20:
           if t==25: return [5.28,2.93,2.04,1.6,1.2,1.0,0.8,0.6,0.5,0.3,0.1],[0.1,0.3,0.5,0.6,0.8,1.0,1.2,1.6,2.0,2.89,5.22]
           elif t==10: return [6.34,3.34,2.27,1.7,1.3,1.0,0.8,0.6,0.4,0.2,0.1],[0.1,0.2,0.4,0.6,0.8,1.0,1.3,1.7,2.22,3.27,6.28]
           elif t==5: return [7.91,3.8,2.49,1.8,1.3,1.0,0.8,0.5,0.3,0.2,0.1],[0.1,0.2,0.3,0.5,0.7,1.0,1.3,1.8,2.46,3.8,7.77]
           elif t==1: return [14.55,5.63,3.33,2.17,1.5,1.0,0.7,0.5,0.3,0.1,0.0],[0.0,0.1,0.2,0.4,0.7,1.0,1.3,2.01,3.07,5.28,13.32]
           elif t==0.5: return [19.18,6.47,3.68,2.32,1.4,1.0,0.7,0.4,0.2,0.1,0.1],[0.0,0.1,0.2,0.4,0.6,1.0,1.2,2.03,3.29,6.0,16.23]
        elif R2==21:
           if t==25: return [5.72,3.0,2.12,1.6,1.3,1.0,0.8,0.6,0.4,0.3,0.1],[0.1,0.3,0.4,0.6,0.8,1.0,1.2,1.6,2.08,3.08,5.78]
           elif t==10: return [7.0,3.5,2.42,1.7,1.3,1.0,0.8,0.6,0.4,0.2,0.1],[0.1,0.2,0.4,0.5,0.7,1.0,1.3,1.7,2.33,3.58,6.99]
           elif t==5: return [9.0,4.0,2.73,1.9,1.4,1.0,0.7,0.5,0.3,0.2,0.1],[0.1,0.2,0.3,0.5,0.7,1.0,1.3,1.8,2.61,4.24,8.72]
           elif t==1: return [15,6.0,3.35,2.17,1.6,1.0,0.7,0.4,0.2,0.1,0.0],[0.0,0.1,0.2,0.5,0.6,1.0,1.4,2.04,3.29,6.37,15.78]
           elif t==0.5: return [20,7.0,3.81,2.4,1.7,1.0,0.7,0.4,0.2,0.1,0.0],[0.0,0.1,0.2,0.4,0.6,1.0,1.5,2.3,3.87,7.77,22.05]
        elif R2==22:
           if t==25: return [5.92,3.13,2.12,1.6,1.3,1.0,0.8,0.6,0.4,0.3,0.1],[0.1,0.3,0.4,0.6,0.8,1.0,1.2,1.6,2.08,3.08,5.78]
           elif t==10: return [7.3,3.68,2.42,1.7,1.3,1.0,0.8,0.6,0.4,0.2,0.1],[0.1,0.2,0.4,0.5,0.7,1.0,1.3,1.7,2.33,3.58,6.99]
           elif t==5: return [9.31,4.36,2.73,1.9,1.4,1.0,0.7,0.5,0.3,0.2,0.1],[0.1,0.2,0.3,0.5,0.7,1.0,1.3,1.8,2.61,4.24,8.72]
           elif t==1: return [16.97,6.35,3.35,2.17,1.6,1.0,0.7,0.4,0.2,0.1,0.0],[0.0,0.1,0.2,0.5,0.6,1.0,1.4,2.04,3.29,6.37,15.78]
           elif t==0.5: return [22.76,7.67,3.81,2.4,1.7,1.0,0.7,0.4,0.2,0.1,0.0],[0.0,0.1,0.2,0.4,0.6,1.0,1.5,2.3,3.87,7.77,22.05]
        return [],[]  



    def draw_preds(self, axes, result, name): 
        xt = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95]
        xn = ['$<5$', '5-15', '15-25', '25-35', '35-45', '45-55','55-65', '65-75', '75-85', '85-95','${>}95$']
        yMax1, yMax2 = [], [] 
        for i,spot in enumerate(['b1','b2','a1','a2']): 
            ax = axes[i] 
            ps = result['preds'][spot] 
            X, Y = ps.X, ps.odds
            ax.plot(X, Y, marker='o') 
            if i in [0,2]: 
                if ps.loc == int(ps.loc): ts = str(int(ps.loc))+'\%' 
                else:                     ts = str(ps.loc)+'\%' 
                ax.set_title('Prediction of Lower '+ts, y = 0.8) 
                E1, E2 = self.get_expected_pred(int(100*result['r2']), ps.loc) 
                if len(E1) == len(X): 
                    ax.plot(X, E1, color = 'orange')  
                    ax.fill_between(X, Y, E1, color = 'orange', alpha=0.25) 
            else: 
                ax.set_title('Prediction of Upper '+ts, y = 0.8) 
                if len(E2) == len(X): 
                    ax.plot(X, E2, color = 'orange') 
                    ax.fill_between(X, Y, E2, color = 'orange', alpha=0.25) 
            for x,(c1,c2) in zip(ps.X, ps.ci): 
                if x == 50: ax.scatter(x, 1, marker='s', color='k', s=50, zorder=10)
                else:       ax.plot([x,x],[c1,c2], color='blue', lw=2) 
            if i < 2: yMax1.append(ax.get_ylim()[1]) 
            else:     yMax2.append(ax.get_ylim()[1]) 
            if i in [2,3]: 
                ax.set_xticks(xt) 
                ax.set_xticklabels(xn,fontsize=8,rotation=50) 
            ax.set_ylabel('Odds Ratio') 
        for ax in axes[0:2]: ax.set_ylim(-1, max(yMax1)) 
        for ax in axes[2::]: ax.set_ylim(0, max(yMax2)*1.2) 





    def draw_curve(self, ax, result, name): 
        X, Y1, Y2 = result['means'] 
        if len(name) < 24: ax.set_title(" ".join(name.split('_')), fontsize=self.fs1-1, x=0.015, y=0.96, ha='left',va='top')
        else:              ax.set_title(" ".join(name.split('_'))[0:20]+'.', fontsize=self.fs1-1, x=0.015, y=0.96, ha='left',va='top')
        ax.plot(X, Y1, color = 'darkorange') 
        ax.scatter(X[1:-1], Y2[1:-1], color='dodgerblue',ec='k') 
        p1, p2 = result['pvals'] 
        e1, e2 = result['effects'] 
        if p1 > 0.05 or e1[0] < 0: ax.scatter(X[0], Y2[0], marker='v', color='dodgerblue', ec='k', s=111) 
        else:                      ax.scatter(X[0], Y2[0], marker='^', color='dodgerblue', ec='k', s=111) 
        if p2 > 0.05 or e2[0] > 0: ax.scatter(X[-1], Y2[-1], marker='^', color='dodgerblue', ec='k', s=111) 
        else:                      ax.scatter(X[-1], Y2[-1], marker='v', color='dodgerblue', ec='k', s=111)        
        ax.set_yticks([]) 
        ax.set_xlabel('Trait Percentile',fontsize=self.fs2) 
        ax.set_ylabel('PRS', fontsize=self.fs2) 
        return



    def finish(self, name): 
        fig_name = self.args.out+'-tailSummary-'+'_'.join(name.split())
        plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.94,wspace=0, hspace=0) 
        if  self.args.plotFormat == 'pdf': plt.savefig(fig_name+'.pdf', dpi=300) 
        else:                              plt.savefig(fig_name+'.png', dpi=300) 
        self.fig.clf() 
            

#############################################################################################
##################################  MULTI TRAIT PLOT     ####################################
#############################################################################################

class MultiPlot:
    def __init__(self, args, progress, results):
        self.args, self.progress, self.results = args, progress, results
        self.filenum = len(self.results.keys()) 
        self.names = [k for k in self.results.keys()]
        self.fs1, self.fs2, self.fs3, self.fps = 34, 24, 15, 23 
        
        if self.filenum   <= 10:   self.fps, self.ds = 19, 200  
        elif self.filenum <= 15: self.fps, self.ds = 18, 166 
        elif self.filenum <= 20: self.fps, self.ds = 16, 133  
        elif self.filenum <= 30: self.fps, self.ds = 14, 100 
        else:                   self.fps, self.ds = 12, 65 
        self.tp = TailPlot(args, progress)  


    def make_multi_curves(self): 
        figNum = 1 
        for gi in range(0,len(self.names),20): 
            group = self.names[gi:gi+20] 
            if len(group) == 1:     self.WD, self.HT, self.rows, self.cols = 10, 8, 1, 1  
            elif len(group) == 2:   self.WD, self.HT, self.rows, self.cols = 17, 8, 1, 2  
            elif len(group) == 3:   self.WD, self.HT, self.rows, self.cols = 24, 8, 1, 3  
            elif len(group) == 4:   self.WD, self.HT, self.rows, self.cols = 17, 15, 2, 2  
            elif len(group) <= 6:   self.WD, self.HT, self.rows, self.cols = 15, 20, 3, 2  
            elif len(group) <= 8:   self.WD, self.HT, self.rows, self.cols = 15, 23, 4, 2  
            elif len(group) <= 9:   self.WD, self.HT, self.rows, self.cols = 16, 20, 3, 3  
            elif len(group) <= 10:  self.WD, self.HT, self.rows, self.cols = 15, 25, 5, 2  
            elif len(group) <= 12:  self.WD, self.HT, self.rows, self.cols = 16, 23, 4, 3  
            elif len(group) <= 15:  self.WD, self.HT, self.rows, self.cols = 16, 25, 5, 3  
            else:                   self.WD, self.HT, self.rows, self.cols = 18, 25, 5, 4  
            
            self.fig = matplotlib.pyplot.gcf()
            self.axes, self.ax_index = [], 0  
            self.fig.set_size_inches(self.WD, self.HT) 
            for i in range(self.rows): 
                for j in range(self.cols):  
                    self.axes.append(plt.subplot2grid((self.rows, self.cols), (i,j),     rowspan=1, colspan=1))
            for i,n in enumerate(group):
                res = self.results[n] 
                self.tp.draw_curve(self.axes[i], res, n)  
            while i + 1 < len(self.axes): 
                i+=1 
                self.axes[i].axis('off') 
            if self.filenum <= 6: fig_name = self.args.out+'-dists'
            else:                 fig_name = self.args.out+'-dists'+str(figNum) 
            figNum+=1 
            if len(group) < 4:     plt.subplots_adjust(left=0.1, bottom=0.15, right=0.95, top=0.94,wspace=0.2, hspace=0.2) 
            elif len(group) < 7:   plt.subplots_adjust(left=0.1, bottom=0.15, right=0.95, top=0.94,wspace=0.2, hspace=0.3) 
            else:                  plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.96,wspace=0.2, hspace=0.3) 
            if  self.args.plotFormat == 'pdf': plt.savefig(fig_name+'.pdf', dpi=300) 
            else:                              plt.savefig(fig_name+'.png', dpi=300) 
            self.fig.clf() 
            


    def make_multi_forests(self): 
        self.fig = matplotlib.pyplot.gcf()
        self.axes, self.ax_index = [], 0  
        hl = int(len(self.names)/4.0)
        self.WD, self.HT, self.rows, self.cols = 10, 3+(hl), 1, 1  
        self.fig.set_size_inches(self.WD, self.HT) 
        self.ax = plt.subplot2grid((self.rows, self.cols), (0,0),     rowspan=1, colspan=1)
        self.yloc = 0 
        for i,n in enumerate(self.names): 
            self.add_forest(n, args.colors[i]) #cmap(i)) 
        self.yloc += 0.5
        self.ax.set_yticks([]) 
        self.ax.set_xticks([-2.5, -2,-1.5,-1,-0.5,0.5,1,1.5,2, 2.5]) 
        self.ax.set_xticklabels([-1*(-2.5*self.rO), -1*(-2+self.rO),-1*(-1.5+self.rO),-1*(-1+self.rO),-1*(-0.5+self.rO),0.5-self.rO,1-self.rO,1.5-self.rO,2-self.rO,2.5-self.rO], fontsize=self.fs3)
        self.ax.set_xlim(-2.55,2.55) 
        self.ax.set_ylim(self.yloc,1) 
        self.ax.set_xlabel('POPout\nEffects', fontsize=self.fs2-2) 
        self.ax.text(-1.5,self.yloc*1.18-0.33,'Lower Tail', ha='center',fontsize=self.fs2-2, fontweight='bold') 
        self.ax.text(1.5,self.yloc*1.18-0.33,'Upper Tail', ha='center',fontsize=self.fs2-2, fontweight='bold') 
        for xl in [-2,-1.5,-1,-0.5,0.5,1,1.5,2]: 
            self.ax.plot([xl,xl],[self.yloc, 1], color='k', linestyle='--', alpha=0.1) 
        fig_name = self.args.out+'-forest'
        plt.subplots_adjust(left=0.05, bottom=0.30, right=0.95, top=0.94,wspace=0.2, hspace=0.2) 
        if  self.args.plotFormat == 'pdf': plt.savefig(fig_name+'.pdf', dpi=300) 
        else:                              plt.savefig(fig_name+'.png', dpi=300) 
        self.fig.clf() 



    def add_forest(self, n, clr): 
        self.rO =  1
        p1, p2 = self.results[n]['pvals'] 
        eL, eL1, eL2  = [-1*(self.rO + x) for x in self.results[n]['effects'][0]] 
        eH, eH1, eH2  = [self.rO + (x * -1) for x in self.results[n]['effects'][1]] 
        QC = self.results[n]['QC'] 
        if len(n) < 16:  self.ax.text(0, self.yloc, " ".join(n.split('_')), ha='center', va='center',fontsize=self.fps) 
        elif len(n) < 18: self.ax.text(0, self.yloc, " ".join(n.split('_')), ha='center', va='center',fontsize=self.fps-1) 
        elif len(n) < 21: self.ax.text(0, self.yloc, " ".join(n.split('_')), ha='center', va='center',fontsize=self.fps-2) 
        else:             self.ax.text(0, self.yloc, " ".join(n.split('_'))[0:16]+'.', ha='center', va='center',fontsize=self.fps-1) 
        self.ax.plot([eL1, eL2], [self.yloc, self.yloc], color = clr, lw = 2, zorder=1) 
        self.ax.plot([eH1, eH2], [self.yloc, self.yloc], color = clr, lw = 2, zorder=1) 
        if not QC: 
            self.ax.scatter(eL, self.yloc, facecolor = 'white', edgecolor = clr, lw = 2, s = self.ds, zorder=2) 
            self.ax.scatter(eH, self.yloc, facecolor='white', edgecolor = clr, lw=2, s = self.ds,zorder=2)  
        else: 
            if self.results[n]['pvals'][0] < 0.05: self.ax.scatter(eL, self.yloc, fc = clr, ec='k', lw=2,s = self.ds) 
            else:                                  self.ax.scatter(eL, self.yloc, color = clr, s = self.ds) 
            if self.results[n]['pvals'][1] < 0.05: self.ax.scatter(eH, self.yloc, fc = clr, ec = 'k', lw=2 ,s = self.ds)  
            else:                                  self.ax.scatter(eH, self.yloc, color = clr, s = self.ds)  
        self.yloc -= 1






#############################################################################################
##################################     POPOUT CLASS      ####################################
#############################################################################################


class PopOdds:
    def __init__(self, loc, tailOdds, fullOdds, SK): 
        self.loc = loc
        self.lowerTailOdds = tailOdds['lt'] 
        self.upperTailOdds = tailOdds['ut'] 
        self.X, self.odds, self.ci = [], [], [] 
        for x in [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]: 
            if x == 'lt':   self.X.append(100*SK) 
            elif x == 'ut': self.X.append(100-100*SK) 
            elif x == 50:   self.X.append(50) 
            elif x < 50:    self.X.append(x+5) 
            else:           self.X.append(x-5) 
            self.odds.append(fullOdds[x][0]) 
            self.ci.append(fullOdds[x][1::])












class PopMunge:
    def __init__(self, args, name, K): 
        self.args, self.name = args, name 
        self.r2 = K['r2'] 
        self.p1, self.p2 = K['pvals'] 
        self.e1, self.e2 = K['effects'][0][0], K['effects'][1][0]
        self.s1, self.s2 = K['se'] 
       
        self.estimates = dd(lambda: {})  
        self.r2 = 0.068
        self.p1, self.p2 = 0.094, 3.8e-39
        self.e1, self.e2 = -0.0761, 0.4963
        self.s1, self.s2 = 0.027565, 0.028855




    def iter_estimate(self, betas, I):
        lower_samples = np.random.normal(loc=self.e1, scale=self.s1, size=I)
        upper_samples = np.random.normal(loc=self.e2, scale=self.s2, size=I)
        lower_samples = np.clip(lower_samples, 0, None)
        upper_samples = np.clip(upper_samples, 0, None)
        for bi,b in enumerate(betas):
            Bk = dd(list)
            for k in range(I):
                p_in_tail = est_prop_in_tail([lower_samples[k], upper_samples[k]], beta=b, r2=self.r2)
                Bk['v1'].append(p_in_tail[0])
                Bk['v2'].append(p_in_tail[1])
                p_in_tail = np.minimum(p_in_tail, 1.0)
                ex = h2_rare_big(p_in_tail, beta=b, rare_maf=1e-5)
                for ek,ev in ex.items(): Bk[ek].append(ev)
                Bk['theta1'].append(p_in_tail[0])
                Bk['theta2'].append(p_in_tail[1])
            for x in ['h2','theta1','theta2']: self.estimates[b][x] = get_quantiles(Bk[x])
            for x in ['m1','m2']: self.estimates[b][x] = np.median(Bk[x])
            for x in ['v1','v2']: self.estimates[b][x] = cond_bool(Bk[x])
        return

    

#############################################################################################
##################################     POPOUT CLASS      ####################################
#############################################################################################



class PopOut:
    def __init__(self, args, progress): 
        self.args, self.progress, self.results, self.names = args, progress, dd(lambda: {}), [] 
        self.SS, self.SK = self.args.tailSize, self.args.tailSize/100.0
        if self.SS > 25: self.PopError('Invalid TailSize (Must be below 25%)')
        if self.SS < 1 and self.SS not in [0.5,0.25,0.1]: self.PopError('Invalid Tailsize - Accepted Fractions are 0.5, 0.25, 0.1') 
        self.MUNGE = False 
        self.R2 = True 

        self.prepare_output_files() 

    def PopError(self, msg):  
        sys.stderr.write('POPError- '+msg+'\n') 
        sys.exit() 
    
    def process(self): 
        if len(self.names) > 1: 
            plot = MultiPlot(self.args, self.progress, self.results) 
            plot.make_multi_curves() 
            plot.make_multi_forests() 

    def prepare_output_files(self): 

        
        self.fh = {'pop': open(self.args.out+'-popout.txt','w')} 
        self.fh['pop'].write('%-30s %8s %8s %8s %9s %9s %10s %10s' % ('---', 'tail','len', 'r2', 'pv1','pv2','se1','se2')) 
        self.fh['pop'].write(' %7s %7s %7s %7s %7s %7s %6s\n' % ('e1', 'e1Lo', 'e1Hi','e2','e2Lo','e2Hi','QC')) 
        if self.R2:
            self.fh['r2'] = open(self.args.out+'-R.txt','w') 
            self.fh['r2'].write('%-30s %8s ' % ('---', 'R-tot')) 
            self.r_locs = [1,5,10,25,75,90,95,99] 
            self.r_format = ' '.join(['%7s' for rl in self.r_locs])+' ' 
            self.fh['r2'].write(self.r_format % tuple(['R'+str(rl) for rl in self.r_locs])) 
            self.fh['r2'].write(self.r_format % tuple(['p'+str(rl) for rl in self.r_locs])) 
            self.fh['r2'].write('\n') 

            #self.fh['r2'].write('%-35s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n' % ('---', 'R-tot','R1', 'R5', 'R10','R25','R75','R90','R95','R99')) 




        if self.args.savePts: self.fh['pts'] = open(self.args.out+'-pts.txt','w') 
        if self.MUNGE: 
            self.fh['munge'] = open(self.args.out+'-h2her.txt','w') 
            self.fh['munge'].write('%-30s %10s %10s %10s %10s %10s %10s\n' % ('---', 'beta','h2', 'h2Lo', 'h2Hi','theta1','theta2')) 
        return 


    
    def write_out(self, MUNGE=True): 
        w, n, r = self.fh['pop'], self.name, self.results[self.name] 
        w.write('%-30s %8s %8d %8.2f' % (n, self.args.tailSize, r['len'], r['r2'])) 
        for p in r['pvals']: 
            if p > 0.001:  w.write(' %9.3f' % p) 
            else:          w.write(' %9.1e' % p) 
        w.write(' %10.5f %10.5f' % (r['se'][0], r['se'][1])) 
        w.write(' %7.2f %7.2f %7.3f' % tuple(r['effects'][0])) 
        w.write(' %7.2f %7.2f %7.3f' % tuple(r['effects'][1])) 
        w.write(' %6s\n' % r['QC']) 
        try: 
            w = self.fh['r2'] 
            w.write('%-30s %8.3f ' % (n, r['R']['tot'][0])) 
            r_vals, p_vals = tuple([str(round(r['R'][k][0],3)) for k in self.r_locs]), [] 
            for k in self.r_locs: 
                if r['R'][k][1] < 0.001: p_vals.append('<0.001') 
                else:                    p_vals.append(str(round(r['R'][k][1],3))) 
            self.fh['r2'].write(self.r_format % r_vals) 
            self.fh['r2'].write(self.r_format % tuple(p_vals)) 
            self.fh['r2'].write('\n') 

        except: pass                 


        try: 
            w = self.fh['pts'] 
            X, Y1, Y2 = r['means'] 
            w.write('%-20s %10s %5s %s\n' % (n, self.args.tailSize, 'X', ','.join([str(x) for x in X]))) 
            w.write('%-20s %10s %5s %s\n' % (n, self.args.tailSize, 'Ye', ','.join([str(x) for x in Y1]))) 
            w.write('%-20s %10s %5s %s\n' % (n, self.args.tailSize, 'Yo', ','.join([str(x) for x in Y2]))) 
        except: pass  
        try: 
            w = self.fh['munge'] 
            for b in self.args.betas: 
                res = r['munge'].estimates[b] 
                w.write('%-30s %10.2f %10.3f %10.3f %10.3f %10.3f %10.3f\n' % (n, b, res['h2'][0], res['h2'][1], res['h2'][2], res['theta1'][0], res['theta2'][0]))
        except: pass  
        




    def read_in(self, f, f_ext): 
        lp = f.readline().split() 
        if len(lp) == 0: self.PopError('Empty POPfile: '+f.name+'\n') 
        if f_ext == 'pop' and len(lp) == 2 and self.args.prs_col == None and self.args.pheno_col == None:  
            self.args.prs_col, self.args.pheno_col = 0, 1  
        else:
            if self.args.prs_col == None and self.args.pheno_col == None: self.PopError('Custom File Format Requires Column Indices (--prs-col, --pheno-col)') 
            if self.args.prs_col == self.args.pheno_col: self.PopError('--prs-col and --pheno-col must be distinct') 
            if self.args.prs_col < 0 or self.args.prs_col >= len(lp): self.PopError('--prs-col must fit within file columns') 
            if self.args.pheno_col < 0 or self.args.pheno_col >= len(lp): self.PopError('--pheno-col must fit within file columns') 
        self.phenoPrs = [] 
        if self.args.no_header: 
            try: prs, pheno = float(lp[self.args.prs_col]), float(lp[self.args.pheno_col]) 
            except ValueError: self.PopError('Invalid Line: '+" ".join(lp)+'\n          Turn off --no-header?')   
            self.phenoPrs.append([pheno,prs]) 
        for i,lp in enumerate(f): 
            line = lp.split()
            try: prs, pheno = float(line[self.args.prs_col]), float(line[self.args.pheno_col]) 
            except ValueError: self.PopError('Invalid Line: '+lp.strip()) 
            self.phenoPrs.append([pheno,prs]) 
        self.sampleLen = len(self.phenoPrs) 
        pStd = np.std([p[1] for p in self.phenoPrs]) 
        self.phenoPrs = [[p[0], p[1], p[1]/pStd] for p in self.phenoPrs] 
        self.phenoPrs.sort(key = lambda X: X[0]) 
        for i,(pheno,prs,prs_std) in enumerate(self.phenoPrs): self.phenoPrs[i].append(min(99.999,round(100*(i/self.sampleLen),3)))
        self.phenoPrs.sort(key = lambda X: X[1]) 
        for i,(pheno,prs, prs_std, pheno_qt) in enumerate(self.phenoPrs): self.phenoPrs[i].append(min(99.999,round(100*(i/self.sampleLen),3)))
        self.phenoPrs.sort() 
        return

    def run_analysis(self, PREDS=True, MUNGE=False):
        self.get_tail_popout() 
        if PREDS: 
            self.progress.begin_analysis(COMMAND='PREDS') 
            self.get_preds() 
            self.progress.step(9) 
        if MUNGE: 
            self.progress.begin_analysis(COMMAND='MUNGE') 
            self.get_munge() 
            self.progress.step(1) 
        if not self.args.no_plots: 
            self.progress.begin_analysis(100*self.SK, COMMAND='PLOT') 
            tp = TailPlot(self.args, self.progress).create_trait_plot(self, self.name) 
        return


    def get_tail_popout(self): 
        self.progress.begin_analysis(100*self.SK, COMMAND='POPOUT') 
        self.results[self.name]['len'] = self.sampleLen 
        self.phenos, self.prs, self.prs_std = [d[0] for d in self.phenoPrs], [d[1] for d in self.phenoPrs], [d[2] for d in self.phenoPrs]
        self.results[self.name]['R'] = {'tot': stats.pearsonr(self.phenos, self.prs)} 
        for job in [self.run_test, self.run_empirical_pv, self.run_qc, self.get_means]: 
            self.progress.step(2) 
            job() 
            self.progress.step(2) 
        return

    def read_and_run(self, i, f, f_ext, name): 
        self.name = name 
        self.names.append(name) 
        self.progress.initialize(f.name, name, i) 
        self.read_in(f, f_ext) 
        self.run_analysis() 
        self.write_out() 

    def get_munge(self): 
        munge = PopMunge(self.args, self.name, self.results[self.name]) 
        munge.iter_estimate(self.args.betas, self.args.iter)
        self.results[self.name]['munge'] = munge 

    def get_preds(self, iters=10): 
        if self.args.predSize == 0: self.preds = [self.args.tailSize] 
        if len(self.preds) == 1: self.preds.append(25) 
        self.preds.sort(reverse=True) 
        a1, a2  = self.preds[0], 100-self.preds[0] 
        b1, b2 =  self.preds[1], 100-self.preds[1] 
        pKey, pLists = {}, dd(lambda: dd(list)) 
        pNames = {'a1': a1, 'a2': a2, 'b1': b1, 'b2': b2} 
        for pheno, prs, prs_std, pheno_qt, prs_qt in self.phenoPrs:
            PK = dd(bool) 
            if pheno_qt < a1: 
                PK['a1'] = True 
                if pheno_qt < b1: PK['b1'] = True 
            elif pheno_qt > a2: 
                PK['a2'] = True 
                if pheno_qt > b2: PK['b2'] = True 
            nt = [int(round(prs_qt,-1))] 
            if prs_qt < b1: nt.append('lt') 
            if prs_qt > b2: nt.append('ut') 
            for k in ['a1','a2','b1','b2']: 
                for nx in nt: pLists[k][nx].append(int(PK[k])) 

        for k in ['a1','a2','b1','b2']: 
            kOdds, notUpper, notLower = {}, dd(float), dd(float) 
            for loc,vals in pLists[k].items(): 
                vSum, vNot, vL = sum(vals), len(vals) - sum(vals), int(len(vals) / 2) 
                for lx,ln in [['ut',notUpper],['lt',notLower]]: 
                    if loc != lx: 
                        ln[0] += vNot 
                        ln[1] += vSum 
                if vNot == 0: 
                    odds = [vSum/0.01, [vSum/0.01 for i in range(iters)]] 
                else: 
                    odds = [vSum / vNot, []] 
                    for i in range(iters): 
                        random.shuffle(vals) 
                        vhalf = vals[0: vL] 
                        try: odds[1].append(sum(vhalf) / (vL - sum(vhalf))) 
                        except ZeroDivisionError: odds[1].append(vL / 0.01) 
                kOdds[loc] = odds 
            
            tailOdds, fullOdds = {}, {} 
            tailOdds['lt'] =  kOdds['lt'][0] / (notLower[1]/notLower[0])
            tailOdds['ut']  = kOdds['ut'][0] / (notUpper[1]/notUpper[0])            
            mScr, mList = kOdds[50][0], [max(0.0001,mx) for mx in kOdds[50][1]] 
            for iLoc, (iScr, iList)  in kOdds.items(): 
                if iLoc == 50: fullOdds[iLoc] = [1, 1, 1] 
                elif iLoc in ['lt','ut']: continue 
                else:  
                    if iScr == 0: fullOdds[iLoc] = [0.01, 0.01, 0.01] 
                    elif mScr == 0: fullOdds[iLoc] = [99.9, 99.9, 99.9] 
                    else: 
                        myVal, myList = iScr/mScr, sorted([iStep/mStep for (iStep, mStep) in zip(iList, mList)])
                        if myVal > 100:    fullOdds[iLoc] = [99.9, min(myList[0],99.9), 99.9] 
                        elif myVal < 0.01: fullOdds[iLoc] = [0.01, 0.01, max(0.01, myList[-2])] 
                        else:              fullOdds[iLoc] = [myVal, max(0.01, myList[1]), min(99.9, myList[-2])]
            pKey[k] = PopOdds(pNames[k], tailOdds, fullOdds, self.SK) 
        self.results[self.name]['preds'] = pKey 






    def get_qt(self, qt): 
        
        if self.SS == 1: return int(qt) 
        elif self.SS > 1: return round(float(int(qt)) / self.SS) * self.SS  
        elif self.SS == 0.5: return round(round(round(qt,1) *2) /2,1)
        elif self.SS == 0.25: return round(round(round(qt,1) *4) /4,1)
        elif self.SS == 0.1:  return round(qt,1)  
        else: print('not supported') 



    def get_means(self): 
        prs_obs, prs_exp = dd(list), dd(list) 
        pObs, pLo = dd(lambda: [[], []]), [1,5,10, 25]  
        if self.SS not in pLo: pLo.append(self.SS) 
        pLo.sort(reverse=True) 
        pHi = [100-p for p in pLo]         
        yInt, beta = self.results[self.name]['params'] 
        for pheno, prs, prs_std, pheno_qt, prs_qt in self.phenoPrs: 
            bin_val = self.get_qt(pheno_qt)
            prs_obs[bin_val].append(prs_std) 
            prs_exp[bin_val].append(yInt + beta*pheno) 
            k = 0 
            if bin_val <= pLo[k]: 
                while k < len(pLo) and bin_val <= pLo[k]: 
                    pObs[pLo[k]][0].append(pheno) 
                    pObs[pLo[k]][1].append(prs_std) 
                    k+=1 
            elif bin_val >= pHi[k]: 
                while k < len(pHi) and bin_val >= pHi[k]: 
                    pObs[pHi[k]][0].append(pheno) 
                    pObs[pHi[k]][1].append(prs_std) 
                    k+=1 
            else: 
                pObs[50][0].append(pheno) 
                pObs[50][1].append(prs_std) 
        for k in pObs.keys(): self.results[self.name]['R'][k] = stats.pearsonr(pObs[k][0], pObs[k][1]) 
        X = [x for x in prs_obs.keys()] 
        prs_obs = [np.mean(prs_obs[x]) for x in X] 
        prs_exp = [np.mean(prs_exp[x]) for x in X] 
        self.results[self.name]['means'] = [X, prs_exp, prs_obs] 
        return 








    def run_test(self):
        
        # Normalize the PRS 
        #ms = np.std(self.prs) 
        #prs = self.prs / ms 
        #self.phenos, self.prs, self.prs_std = [d[0] for d in self.phenoPrs], [d[1] for d in self.phenoPrs], [d[2] for d in self.phenoPrs]
        
        prs = [d[2] for d in self.phenoPrs]
        # Fit the regression model
        X = sm.add_constant(self.phenos)  # Add a constant to the predictor variable (for the intercept)
        self.model = sm.OLS(prs, X).fit()
        self.results[self.name]['params'] = self.model.params 
        R, pv = stats.pearsonr(prs, self.phenos) 
        self.results[self.name]['r2'] = R*R 
        # Test upper tail "regression to mean"
        T_upper = np.percentile(self.phenos, 100 - self.SK * 100)
        upper_tail_indices = np.where(self.phenos > T_upper)[0]
        residuals_upper = self.model.resid[upper_tail_indices]
        upper_se =  np.std(residuals_upper) / np.sqrt(len(residuals_upper))
        t_stat_upper, p_value_upper = stats.ttest_1samp(residuals_upper, 0)
        ci_upper = stats.t.interval(0.95, len(residuals_upper) - 1, loc=np.mean(residuals_upper), scale=np.std(residuals_upper) / np.sqrt(len(residuals_upper)))
        # Test lower tail "regression to mean"
        T_lower = np.percentile(self.phenos, self.SK * 100)
        lower_tail_indices = np.where(self.phenos < T_lower)[0]
        residuals_lower = self.model.resid[lower_tail_indices]
        lower_se =  np.std(residuals_lower) / np.sqrt(len(residuals_lower))
        t_stat_lower, p_value_lower = stats.ttest_1samp(residuals_lower, 0)
        ci_lower = stats.t.interval(0.95, len(residuals_lower) - 1, loc=np.mean(residuals_lower), scale=np.std(residuals_lower) / np.sqrt(len(residuals_lower)))
        # Standard deviation of prs
        # Organizing the result
        self.results[self.name]['pvals'] = [p_value_lower, p_value_upper] 
        self.results[self.name]['se'] = [lower_se, upper_se] 
        self.results[self.name]['effects'] = [[np.mean(residuals_lower), ci_lower[0], ci_lower[1]]]
        self.results[self.name]['effects'].append([np.mean(residuals_upper), ci_upper[0], ci_upper[1]]) 
        return











    def run_empirical_pv(self, n_q=100):

        resids = dd(list) 
        yInt, beta = self.model.params  
        for pheno, prs, prs_std, pheno_qt, prs_qt in self.phenoPrs: 
            pred = yInt + beta*pheno 
            bin_val = self.get_qt(pheno_qt) 
            resids[bin_val].append(pred-prs) 
            #resids[int(pheno_qt)].append(pred-prs) 
        emps, emp_ranks = [], [] 
        for k,R in resids.items(): 
            t_stat, pval = stats.ttest_1samp(R,0) 
            emps.append([t_stat,pval,k]) 
        for i,(ts,pv,loc) in enumerate(sorted(emps)): 
            emp_ranks.append([loc,i,ts,pv]) 
        emp_ranks.sort()  
        p = (1+emp_ranks[0][1])/len(emp_ranks), (len(emp_ranks) - emp_ranks[-1][1])/len(emp_ranks)
        self.results[self.name]['empirical'] = p 




    def run_qc(self,stepsize=1): 
        # se 
        trunc_data = [d for d in self.phenoPrs if d[3] > 10 and d[3] < 90] 
        X = sm.add_constant([t[0] for t in trunc_data])  # Add a constant to the predictor variable (for the intercept)
        qc_model = sm.OLS([t[1] for t in trunc_data], X).fit()
        yInt, beta = qc_model.params  
        resids = dd(list) 
        for pheno, prs, prs_std, pheno_qt, prs_qt in trunc_data: 
            pred = yInt + beta*pheno 
            bin_val = self.get_qt(pheno_qt) 
            resids[bin_val].append(pred-prs) 
            #resids[int(pheno_qt)].append(pred-prs) 
        middle_pvs = [] 
        for k,R in resids.items(): 
            t_stat, pval = stats.ttest_1samp(R,0) 
            middle_pvs.append(pval) 
        qc_val = min(middle_pvs) * len(middle_pvs) 
        self.results[self.name]['qc_val'] = qc_val 
        self.results[self.name]['QC'] = qc_val > 0.05 
        return



def load_config(filenames, args, F): 
    cmap = plt.cm.get_cmap('hsv', len(filenames)) 
    if F: 
        N1, N2, N3, NK, CK = {}, {}, {}, {}, {} 
        names, colors = [], [] 
        header = F.readline() 
        if len(header.split(',')) >  1:   header, COMMA = [h.lower() for h in header.strip().split(',')], True  
        else:                             header, COMMA = [h.lower() for h in header.split()], False 
        for i,n in enumerate(filenames): 
            N1[n] = n 
            N2[n.split('/')[-1]] = n
            N3[n.split('/')[-1].split('.pop')[0]] = n 
        if 'name' not in header and 'color' not in header:  
            sys.stderr.write('POPoutError- '+'Invalid Config File - name or color header required'+'\n') 
            sys.exit() 
        for line in F: 
            if COMMA: line = line.strip().split(',') 
            else:     line = line.split() 
            if len(header)   == len(line):  my_header = header[1::] 
            elif len(header) == len(line) - 1: my_header = header 
            else:                              sys.stderr.write('POPoutError- '+'Invalid Config File - Inconsistent Column Numbers'+'\n') 
            K = {h: line[j+1] for j,h in enumerate(my_header)} 
            if line[0] in N1:   nd = N1[line[0]] 
            elif line[0] in N2: nd = N2[line[0]] 
            elif line[0] in N3: nd = N3[line[0]] 
            else:               continue 
            if 'name' in K: NK[nd] = K['name'] 
            if 'color' in K: CK[nd] = K['color'] 
        for i,n in enumerate(filenames):  
            if n in NK: names.append(NK[n])  
            else:       names.append(n.split('/')[-1].split('.pop')[0]) 
        if len(CK) > 0: 
            for i,n in enumerate(filenames): 
                if n in CK: colors.append(CK[n]) 
                else:       colors.append('gray') 
        else:   colors = [cmap(j) for j in range(len(filenames))] 
        return names, colors 
    else: 
        if len(args.names) == len(filenames): names = args.names 
        else:                                 names = [p.name.split('/')[-1].split('.')[0] for p in filenames] 
        if len(args.colors) == len(filenames): colors = args.colors 
        else:                                  colors = [cmap(j) for j in range(len(filenames))] 
        return names, colors 




#def run_script(args, parser, command_line):
def main(args, parser, command_line):

    
    if len(args.popfiles) > 0: 
        progress = PopProgress(args,command_line) 
        popOut = PopOut(args, progress) 
        args.names, args.colors = load_config([p.name for p in args.popfiles], args, args.config) 
        for i,f in enumerate(args.popfiles): popOut.read_and_run(i,f, f.name.split('.')[-1], args.names[i])
        popOut.process() 
        progress.complete_analysis()  
    elif args.makePOPfile: 
        progress = PopProgress(args,command_line, MODE='Make Popfile') 
        K = {} 
        for i,f in enumerate(args.makePOPfile): 
            for j,line in enumerate(f): 
                line = line.split()
                if len(line) == 1: line = line[1].split(',')
                try: xId, xVal = line[0], float(line[-1])  
                except ValueError: continue 
                if i == 0: K[xId] = [xVal]  
                elif xId in K: K[xId].append(xVal) 
                else: continue
        w = open(args.out+'.pop','w') 
        w.write('%-20s %20s\n' % ('PRS','Phenotype')) 
        for V in K.values(): 
            if len(V) == 2: w.write('%-20s %20s\n' % (V[0], V[1]))  
        w.close() 
        progress.complete_analysis() 
    else: 
        parser.print_help() 
    sys.exit() 



if __name__ == '__main__':
    import argparse, sys
    usage = "usage: ./%prog [options] data_file"
    parser = argparse.ArgumentParser()
    parser.allow_abbrev=True
    

    # PRS ON PHENOTYPE 
    parser.add_argument('popfiles',type=argparse.FileType('r'),nargs='*') 
    parser.add_argument('--config',type=argparse.FileType('r')) 
    parser.add_argument("--names", nargs='+',default=[],type=str,metavar='',help="Trait Name(s)")
    parser.add_argument("--colors", nargs='+',default=[],type=str,metavar='',help="Trait Colors (For Forest Plot)")
    parser.add_argument("--tailSize", type=float,default=1.0,metavar='',help="Tail Size (default 1 Percent)")
    parser.add_argument("--predSize", type=float,default=0,metavar='',help="Tail Size (default 1 Percent)")
    parser.add_argument('--makePOPfile',type=argparse.FileType('r'), nargs = 2, metavar='', help='Make POP File (PRS on Phenotype)') 
    parser.add_argument("-o","--out", type=str,default='out',metavar='',help="Output Prefix")
    parser.add_argument("--plotFormat", type=str,default='pdf',metavar='',help="pdf or png")
    parser.add_argument("--savePts", action='store_true', default=False,help="Save Plot Distribution")  
    parser.add_argument("--silent", action='store_true', default=False,help="Suppress Output Stream") 
    parser.add_argument("--no-header", action='store_true', default=False,help="POP file has no header") 
    parser.add_argument("--no-plots", action='store_true', default=False,help="Suppress Plotting") 
    parser.add_argument("--prs-col", default=None,type=int,metavar='',help="PRS Column")
    parser.add_argument("--pheno-col", default=None,type=int,metavar='',help="Phenotype Column")
    #parser.add_argument("--betas", nargs='+',default=[3],type=float,metavar='',help="Trait H2")
    #parser.add_argument("--iter", type=int,default=10,metavar='',help="Iterations")   
    args = parser.parse_args() 
    main(args, parser, ' '.join(sys.argv)) 






