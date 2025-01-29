#!/usr/bin/python3

import matplotlib 
import matplotlib.pyplot as plt
import sys
from collections import defaultdict as dd
import numpy as np
import scipy.stats as stats
import statsmodels.api as sm
import random

#import matplotlib
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes      
matplotlib.rcParams['xtick.labelsize'] = 24                                                                                                                                                                                 
matplotlib.rcParams['ytick.labelsize'] = 24                                                                                                                                                                                 
#matplotlib.rcParams['xtick.major.size'] = 18                                                                                                                                                                               
#matplotlib.rcParams['ytick.major.size'] = 18                                                                                                                                                                               
matplotlib.rcParams['axes.linewidth'] = 1.85  


plt.rc('text', usetex=True)
matplotlib.use("Cairo")
import warnings
warnings.filterwarnings("ignore")




def load_p_lists(PK, p_mean, FULL=True): 
    if FULL: aS,zS = 100, 0 
    else:    aS,zS = 25, 75
    a, z = 0, 100 
    aList, zList, aRes, zRes = [], [], [], [] 
    while a < aS: 
        a+=0.1
        a = round(a,1) 
        aList.extend(PK[a]) 
        aMean = np.mean(aList) 
        T,pv = stats.ttest_1samp(aList, p_mean) 
        if aMean < p_mean: V=True 
        else:              V=False 
        #pd = pdiff(aMean, p_mean) 
        #if pd < 0: V=True 
        #else:      V=False  
        #aRes.append([a, aMean-p_mean, aMean, V, T, pv]) 
        aRes.append([a, aMean, T, pv]) 
    while z > zS: 
        z-=0.1
        z=round(z,1) 
        zList.extend(PK[z]) 
        zMean = np.mean(zList) 
        T,pv = stats.ttest_1samp(zList, p_mean) 
        if zMean < p_mean: V=False 
        else:              V=True 
        zRes.append([z,zMean,T,pv]) 
        #pd = pdiff(zMean, p_mean) 
        #if pd > 0: V=True 
        #else:      V=False  
        #zRes.append([z,pd,zMean,V,T,pv]) 
    return aRes, zRes



def pdiff(a,p):
    
    a = a*10000 
    p = p*10000

    a*=-1 

    return a-p 
    

    print(a,p) 

    sys.exit() 


    pd = round(100*(abs(a-p) / abs((a+p)/2)),2)
    if a < p: return -1*pd 
    else:     return pd 




class MyFigure:
    def __init__(self,fn=6): 
        self.setup(fn) 

    # SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP #
    def setup(self,row_len): 
        
        self.ax_index = 0 
        self.fig, self.axes = matplotlib.pyplot.gcf(), [] 
        self.rows, self.cols, self.WD, self.HT = 4, 4, 50, 50 
        for i in range(self.rows):
            for j in range(self.cols): 
                self.axes.append(plt.subplot2grid((self.rows,self.cols), (i,j), rowspan = 1, colspan=1)) 
        self.fig.set_size_inches(self.WD, self.HT) 
        return         




    def finish(self): 

        letters = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v'] 
        letters.extend(['w','x','y','z','z1','z2']) 

        j = 0 
        for i,ax in enumerate(self.axes): 
            #if i in [4,19]: continue
            lt = letters[j] 
            j+=1 
            (a,b),(c,d) = ax.get_xlim(), ax.get_ylim() 
            xs, ys = (b-a)/30.0, (d-c)/50.0 
            a -= xs 
            d += ys 

            ax.text(a,d,lt,fontsize=35,fontweight='bold',ha='left',va='bottom',clip_on=False) 
            
            if j == len(letters): break 


        
        plt.subplots_adjust(left=0.03, bottom=0.05, right=0.98, top=0.93,wspace=0.15, hspace=0.70) 
        plt.savefig('supFigure2.pdf',dpi=500) 
        plt.savefig('tmpFigure.pdf',dpi=500)  
        #plt.savefig('checkplot_'+mn+'.pdf',dpi=500) 
        #plt.savefig(self.fprefix+'.png',dpi=400)  
        #plt.savefig(self.fprefix+'.pdf',dpi=600)  
        plt.clf() 
        return  


    def add_label(self, ax,pz,c1,c2): 
        x1 = 41
        y1, y2, y3 = 2.15, 1.95,1.60
        ax.text(x1+2,y1,'Case Defns',va='bottom',fontsize=35) 
        ax.plot([x1-3,x1+31],[y1,y1],color='k') 
        ax.text(x1+5,y2-0.1,'Upper '+pz,fontsize=32) 
        ax.plot([x1-1,x1+3],[y2,y2],color=c2) 
        ax.scatter([x1-1,x1+3],[y2,y2],color=c2) 
        ax.text(x1+5,y3-0.08,'Lower '+pz,fontsize=32) 
        ax.plot([x1-1,x1+3],[y3,y3],color=c1) 
        ax.scatter([x1-1,x1+3],[y3,y3],color=c1) 




    def create(self, MD, PK, MEANS, FILE_NAMES): 
        
        trait_ids = [k for k in PK.keys()]
        
        for ti in trait_ids:  
            try: Tnn = NKK[ti] 
            except: Tnn = NK[ti] 
            
            #self.p_analyze(PK[ti], MEANS[ti], FILE_NAMES[ti], self.axes[self.ax_index]) 
            
            #continue 
            #sys.exit() 
            axes = self.axes[self.ax_index: self.ax_index+4] 
            self.md_plot(ti, Tnn, MD, axes) 
            self.ax_index += 4 

        #self.p_analyze(PK, MEANS) 
        #self.md_plot(MD) 
    
    def p_analyze(self, PK, x_means, fn, ax): 
        
        p_mean, t_mean = x_means 
        aRes, zRes = load_p_lists(PK, p_mean) # self.p_mean)  
        
        maxA,maxZ = max([abs(r[1]) for r in aRes]), max([abs(r[1]) for r in zRes]) 
        maxD = max([abs(r[1]) for r in aRes] + [abs(r[1]) for r in zRes])
        aX,aY = [], [] 
        zX,zY = [], [] 


        w = sys.stdout 

        d1 = [fn.split('.')[0], ".".join(fn.split('.')[1::])] 


        aData, zData = [], []  
        for r in aRes: 
            if r[3] < 0.001: pv = "%7.2e" % r[3] 
            else:            pv = str(round(r[3],4)) 
            rk, scr, tv = str(r[0]), str(r[1]), str(round(r[2],4))  
            aData.append(",".join([rk, scr, tv, pv])) 
            aX.append(r[0]) 
            aY.append((r[1]-p_mean)/maxA) 

        for r in zRes: 
            if r[3] < 0.001: pv = "%7.2e" % r[3] 
            else:            pv = str(round(r[3],4)) 
            rk, scr, tv = str(r[0]), str(r[1]), str(round(r[2],4))  
            zData.append(",".join([rk, scr, tv, pv])) 
            zX.append(r[0]) 
            zY.append((r[1]-p_mean)/maxZ) 


        w.write('%-9s %9s %15s %25s %25s %10s' % (d1[0], d1[1], p_mean, maxA, maxZ, 'LOWER:')) 
        w.write(' %s\n' % "|".join(aData))  


        w.write('%-9s %9s %15s %25s %25s %10s' % (d1[0], d1[1], p_mean, maxA, maxZ, 'UPPER:')) 
        w.write(' %s\n' % "|".join(zData))  

        return 

        ax.plot(aX,aY) 
        ax.plot(zX,zY) 
        
        ax.plot([0,100],[0,0],color='k',linestyle='--') 



        #for x,pd,X,T,pv in aRes: 
        #    ax.scatter(x,pd) 



        #[10.0, True, -58.56, -8.321855055870515, 1.329358832730926e-16]]

        return





        print(aRes) 




        sys.exit() 



    def md_plot(self, k, Tnn, MD, axes = 'NA'): 
        X = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95]
        Xl = ['0-5', '5-15', '15-25', '25-35', '35-45', '45-55', '55-65', '65-75', '75-85', '85-95', '95-100']
        X = [5, 10, 15, 20, 25, 30, 35, 40, 45, 52.5, 60, 65, 70, 75, 80, 85, 90, 95, 100] 
        Xl = ['0-5','5-10','10-15','15-20','20-25','25-30','30-35','35-40','40-45','45-55','55-60','60-65','65-70','70-75','75-80','80-85','85-90','90-95','95-100']
        my_colors = ['blue','orange'] 
        for i,o_grp in enumerate([['Lower 10','Upper 10'],['Lower 5','Upper 5'],['Lower 1','Upper 1'],['Lower 0.5','Upper 0.5']]): 
            ax = axes[i] 
            if i == 1: ax.set_title(Tnn, fontsize=88,y=1.20,x=1.1,ha='center') 
            self.add_label(ax,o_grp[0].split()[-1]+'%',my_colors[0],my_colors[1]) 
            for j, (clr,opt) in enumerate(zip(my_colors,o_grp)): 
                Y,Yp = [MD[k][opt][y] for y in X],[MD[k][opt][y][0] for y in X] 
                ax.plot(X[0:9],Yp[0:9],color=clr) 
                ax.plot(X[10::],Yp[10::],color=clr)
                ax.plot(X[8:11],Yp[8:11],color=clr,linestyle='--')
                for xi,(x,[y,yL,yH]) in enumerate(zip(X,Y)): 
                    if x == 52.5: ax.scatter(x,y,marker='s',color='black',s=250, zorder=20)
                    else:       
                        ax.scatter(x,y,color=clr) 
                        ax.plot([x,x],[yL,yH],zorder=20,color=clr) 
                ax.set_xticks(X) 
                ax.set_xticklabels(Xl,rotation=60,fontsize=20) 
                ax.set_xlabel('PRS Quantiles',fontsize=30) 
                ax.set_ylabel('Log(OR)',fontsize=30)  
                ax.set_yticks([-3,-2,-1,0,1,2,3],fontsize=20) 
                ax.set_xlim(3,102) 
                ax.set_ylim(-2.6,2.6) 
        return













































































def zrnd(x):
    if x < 5: return 5
    elif x < 15: return 10
    elif x < 25: return 20 
    elif x < 35: return 30 
    elif x < 45: return 40 
    elif x < 55: return 50
    elif x < 65: return 60 
    elif x < 75: return 70 
    elif x < 85: return 80 
    elif x < 95: return 90 
    else:        return 95 

def mrnd(x):
    if   x < 5: return 5 
    elif x < 10: return 10 
    elif x < 15: return 15 
    elif x < 20: return 20 
    elif x < 25: return 25 
    elif x < 30: return 30 
    elif x < 35: return 35 
    elif x < 40: return 40 
    elif x < 45: return 45 
    #elif x < 50: return 50 
    elif x < 55: return 52.5 
    elif x < 60: return 60 
    elif x < 65: return 65 
    elif x < 70: return 70 
    elif x < 75: return 75 
    elif x < 80: return 80 
    elif x < 85: return 85 
    elif x < 90: return 90 
    elif x < 95: return 95 
    else:        return 100


LKK = {} 
PKK = {} 

LKK[0.5] = ['Lower 0.5','Lower 1','Lower 5','Lower 10']  
LKK[1] = ['Lower 1','Lower 5','Lower 10']  
LKK[5] = ['Lower 5','Lower 10']  
LKK[10] = ['Lower 10']  
#LKK[25] = ['Lower 25','Lower 50']  
#LKK[50] = ['Lower 50']  

#PKK[50] = ['Upper 50']  
#PKK[75] = ['Upper 50','Upper 25']  
PKK[90] = ['Upper 10']  
PKK[95] = ['Upper 10','Upper 5']  
PKK[99] = ['Upper 10','Upper 5','Upper 1']  
PKK[99.5] = ['Upper 10','Upper 5','Upper 1','Upper 0.5']

ALL_KKK = LKK[0.5] + PKK[99.5] 







def rexp(y): 
    if y == 'ALL': return 1.0 
    else:  
        try: 
            z = int(y.split()[-1]) 
            return z/100.0 
        except ValueError: return 0.001  



def recalc(K): 
    
    #div = pt / 100.0 
    T = dd(float) 
    N = dd(lambda: dd(int))     
    
    ODDS = dd(lambda: dd(lambda: [0.0,0.0])) 

    for x in K:
        #T[x] = float(len(K[x]))
        for y in K[x]: 
            if y >= 10 and y < 90: mz = [] 
            else: 
                if y < 0.5: mz = LKK[0.5] 
                elif y < 1: mz = LKK[1] 
                elif y < 5: mz = LKK[5] 
                elif y < 10: mz = LKK[10] 
                elif y > 99.5: mz = PKK[99.5] 
                elif y >= 99: mz = PKK[99] 
                elif y >= 95: mz = PKK[95] 
                elif y >= 90: mz = PKK[90] 
            nz = [z for z in ALL_KKK if z not in mz] 
            for z in mz:
                T[z] += 1  
                N[x][z] += 1 
            
            for z in mz: ODDS[z][x][1] += 1 
            for z in nz: ODDS[z][x][0] += 1 
    
 
    OR = dd(lambda: dd(float))
    RG = dd(lambda: {}) 
    PG = dd(lambda: {}) 
    
    for z in ODDS: 
        Y_50 = [0 for x in range(int(ODDS[z][52.5][0]))] + [1 for x in range(int(ODDS[z][52.5][1]))]
        X_50 = [0 for x in Y_50] 
        for k in sorted(ODDS[z]): 
            
            my_Y = [0 for x in range(int(ODDS[z][k][0]))] + [1 for x in range(int(ODDS[z][k][1]))]
            my_X = [1 for x in my_Y] 
            Y = Y_50 + my_Y 
            X = X_50 + my_X 
            d = {'X': X_50 + my_X, 'Y': Y_50 + my_Y, 'INT': [1 for x in X_50 + my_X]} 
            df = pd.DataFrame(data = d) 
            exog = sm.add_constant(df['X']) 
            glm = sm.GLM(df['Y'], exog, family = sm.families.Binomial()).fit() 
            ci = glm.conf_int(alpha=0.05) 
            PG[z][k] = [glm.params[1], ci[0][1], ci[1][1]] 
    return PG







def run_script(tfiles, pt = 5): 
    MD, MR, AD, MF = {}, {}, {} , {} 
    for tf in tfiles: 
        tn = tf.split('/')[-1].split('.')[0] 
        

        MF[int(tn)] = tf.split('/')[-1]  


        f = open(tf) 
        f.readline() 
        K, P = dd(list), dd(list) 
        for i,line in enumerate(f): 
            line = line.split() 
            p = float(line[1]) 
            if i == 0: 
                p_mean, t_mean = float(line[-2]), float(line[-1]) 
                p_min, p_max = p, p  
            scr = 100*stats.norm.cdf(float(line[2])) 
            
            #print(line) 
            #print(scr) 


            v = round(scr,1) 
            if v == 0: v = 0.1 
            if v == 100: v = 99.9
            if p < p_min: p_min = p 
            if p > p_max: p_max = p 
            P[v].append(p) 
            K[mrnd(int(line[4]))].append(scr) 
        

        MD[int(tn)] = recalc(K)  
        MR[int(tn)] = P 
        AD[int(tn)] = [p_mean, t_mean] # p_min, p_max, t_mean] 

    mf = MyFigure(len(tfiles)) 
    mf.create(MD, MR, AD, MF)
    mf.finish() 




'''
if __name__ == '__main__':
    from optparse import OptionParser
    usage = "usage: ./%prog [options] data_file"
    parser = OptionParser(usage=usage)
    parser.add_option("--pt", default = 5, type=int, help="Output Filename Prefix")
    (options, args) = parser.parse_args()
    run_script(args, options.pt) 
''' 











def sibError(msg): 
    sys.stderr.write('SibError: '+msg+'\n') 
    sys.exit() 

def rage_color_bar(hmap, ax, minv=0,maxv=99,LOW=False,MID=False):                                                                                                                                                                               
    CBAR = plt.cm.ScalarMappable(cmap=hmap, norm=plt.Normalize(vmin = minv, vmax=maxv))                                                                                                                                     
    CBAR._A = []                                                                                                                                                                                                            
    if LOW: axins = inset_axes(ax, width = "200%", height="35%", loc = 'upper left',bbox_to_anchor=(0.66,-1.0,0.4,0.3),bbox_transform=ax.transAxes,borderpad=0,)                                                                    
    elif MID:   axins = inset_axes(ax, width = "200%", height="35%", loc = 'upper left',bbox_to_anchor=(0.09,-0.80,0.4,0.3),bbox_transform=ax.transAxes,borderpad=0,)                                                                    
    else:   axins = inset_axes(ax, width = "200%", height="35%", loc = 'upper left',bbox_to_anchor=(0.10,-0.6,0.4,0.3),bbox_transform=ax.transAxes,borderpad=0,)                                                                    
    plt.colorbar(CBAR, cax = axins, orientation='horizontal', ticks=[])                                                                                                                                                     
    yp = -22                                                                     
    if LOW: 
        ax.text(85,yp-42,'Polygenic Heritability ($h_p^2$)',clip_on=False, fontsize=25)                                                                                                                                          
        ax.text(65,yp-42,'0\%',clip_on=False, fontsize=22)                                                                                                                                                                         
        ax.text(140,yp-41,'100\%',clip_on=False, fontsize=22)     
    elif MID: 
        ax.text(35,yp-23,'Polygenic Heritability ($h_p^2$)',clip_on=False, fontsize=21)                                                                                                                                          
        ax.text(10,yp-26,'0\%',clip_on=False, fontsize=20)                                                                                                                                                                         
        ax.text(85,yp-26,'100\%',clip_on=False, fontsize=20)     
    else: 
        ax.text(28,yp-2,'Polygenic Heritability ($h_p^2$)',clip_on=False, fontsize=25)                                                                                                                                          
        ax.text(10,yp,'0\%',clip_on=False, fontsize=22)                                                                                                                                                                         
        ax.text(83,yp,'100\%',clip_on=False, fontsize=22)     



class SibTable:
    def __init__(self, args, ax,  clr='white', PETT=False): 
        self.args, self.ax = args, ax 
        self.rows, self.cols = [], [] 

    def get_loc(self,X,Y):
        x1,x2 = X[0]/100.0 , X[1]/100.0
        y1,y2 = Y[0]/100.0 , Y[1]/100.0
        return (x1,y1,x2-x1,y2-y1)

    def add_row(self,row_data,X=None,Y=None,COLORS=[],WIDTHS=[],FS=13,ALPHA=0,TITLE=False, CL = 'center',CLEAR=False): 
        if X == None: X = self.X_SPAN
        if Y == None: Y = self.Y_SPAN
        cL,rL,rD,xL = CL,None,[row_data],len(row_data)  
        bl = self.get_loc(X,Y) 
        while len(WIDTHS) < len(row_data): WIDTHS.append(10) 
        while len(COLORS) < len(row_data): COLORS.append('white') 
        if CL != 'center': row = self.ax.table(cellText=rD,rowLabels=rL,cellColours = [COLORS[0:xL]],colWidths=WIDTHS[0:xL],bbox=bl,loc = cL, cellLoc=cL, alpha=ALPHA, clip_on=False) 
        else:              row = self.ax.table(cellText=rD,rowLabels=rL,cellColours = [COLORS[0:xL]],colWidths=WIDTHS[0:xL],bbox=bl,cellLoc=cL, alpha=ALPHA, clip_on=False) 
        row.auto_set_font_size(False)
        row.set_fontsize(FS)
        table_props = row.properties()
        if ALPHA > 0: 
            for cell in row._cells: row._cells[cell].set_alpha(ALPHA)
        self.rows.append(row) 


    def make_res(self, tt, mode, flen): 
        
        if   mode == 1:   fs0,fs1,fs2,fs3 = 16,18,16,14 
        elif mode == 2:   fs0,fs1,fs2,fs3 = 17,14,13,11 
        elif mode == 3:   fs0,fs1,fs2,fs3 = 15,13,11,9 
        else:             fs0,fs1,fs2,fs3 = 16,14,12,10 
        
        c1,c2 = 'lightgrey','whitesmoke'
        SW = [18,14,14]
        th = 9 
        self.add_row(['Lower 1\%','Tail Result','Upper 1\%'],COLORS=[c1,c1,c1],X=(0,100),Y=(88,100), FS = fs0, WIDTHS=[44,12,44], TITLE=True) 
        
        MM = [['Summary'], ['De Novo'], ['Mendelian'], ['Optimal\nRates']]
        ST = [['Samples','Index Mean','Dist-P-val'],['Expected Mean($s_2$)','Observed','P-val'],['Exp Count($s_2 \in$ 1\%)','Observed','P-val'],['De Novo','Mendelian','Polygenic']]
        yL,yS = 88, 22
        for M,S in zip(MM,ST): 
            self.add_row(M,X=(44,56),Y=(yL-yS,yL), WIDTHS=[12], COLORS=[c2],FS=fs1,TITLE=True) 
            self.add_row(S,X=(0,44),Y=(yL-th,yL), WIDTHS=SW, COLORS=[c2,c2,c2],FS=fs2-1,TITLE=True) 
            self.add_row(S,X=(56,100),Y=(yL-th,yL), WIDTHS=SW, COLORS=[c2,c2,c2],FS=fs2-1,TITLE=True) 
            yL -= yS
        rd = [tt['0-0'],tt['99-99']] 
        D1 = [[r.size,round(r.index_avg,2), r.d_str] for r in rd] 
        D2 = [[round(r.n_exp,3), round(r.n_obs,3), r.n_str] for r in rd] 
        D3 = [[round(r.m_exp,2), int(r.m_obs), r.m_str] for r in rd] 
        D4 = [[r.rates['novo'],r.rates['mend'],r.rates['poly']] for r in rd] 
        yL,yS = 88, 22
        for D in [D1,D2,D3,D4]: 
            self.add_row(D[0],X=(0,44),Y=(yL-yS,yL-th), WIDTHS=SW, FS=fs2+1,TITLE=False) 
            self.add_row(D[1],X=(56,100),Y=(yL-yS,yL-th), WIDTHS=SW, FS=fs2+1,TITLE=False) 
            yL-=yS 
        ek = 'darkslategray'
        self.ax.plot([0,100],[101,101],color=ek,linewidth=3,clip_on=False)
        self.ax.plot([0,100],[0,0],color=ek,linewidth=3,clip_on=False)
        self.ax.plot([0,0],[0,101],color=ek,linewidth=3,clip_on=False)
        self.ax.plot([100,100],[0,101],color=ek,linewidth=3,clip_on=False)
        self.ax.set_ylim(0,100) 
        self.ax.set_xlim(0,100) 
        self.ax.axis('off') 
        return self 
    


    def make_h2(self, h2, mode, flen):


        if   mode == 1:   fs0,fs1,fs2,fs3 = 16,15,13.5,12 
        elif mode == 2:   fs0,fs1,fs2,fs3 = 17,14,13,11 
        elif mode == 3:   fs0,fs1,fs2,fs3 = 14,13,11,9 
        else:             fs0,fs1,fs2,fs3 = 16.5,13,11.5,10 
        c1,c2 = 'lightgrey','whitesmoke'
        self.add_row(['Conditional\nHeritability'],COLORS=[c1,c2,c1],X=(0,100),Y=(90,100), FS = fs0, WIDTHS=[100], TITLE=True) 
        self.add_row(['Range\nCovered\n[Samples]','Iterative\nEstimate\n'+'(95\% CI)'],COLORS=[c2,c2],X=(0,100),Y=(77,90), FS = fs1, WIDTHS=[50,50], TITLE=True)  
        h_names, h_locs =    ['All','Bod','Mid','LoH','HiH','LoT','HiT'], ['0-100','5-95','35-65','5-40','60-95','0-4','96-100']
        h_smart = ['Full','Dist Body','Middle','Lower','Upper','Low Tail','Top Tail'] 
        yL,yS = 77, 11
        for i,(k,l,s) in enumerate(zip(h_names, h_locs, h_smart)):
            if i % 2 != 0: myclr = c2 
            else:          myclr = 'white' 
            hd = h2.CI[k] 
            hz = str(h2.PW[k][0])
            h_str = str(hd[0])+'\n('+str(hd[1])+','+str(hd[2])+')'
            self.add_row([s+':\n'+l+'\n['+hz+']',h_str],COLORS=[myclr,myclr],X=(0,100),Y=(yL-yS,yL), FS = fs2, WIDTHS=[50,50], TITLE=True) 
            yL -= yS 
        self.ax.axis('off') 
        return




#############################################################################################
##################################  PROGRESS STREAM  ########################################
#############################################################################################


class PopProgress:
    def __init__(self, args, command_line): 
        self.args = args
        if args.silent: self.ACTIVE = False 
        else:           self.ACTIVE = True 
        #self.out1  = open(self.args.out+'.progress.log','w') 
        self.out2  = sys.stderr 
        self.space = '' 
        self.show('\nPopOut Begins:  '+command_line+'\n')
        self.show('         Mode:  '+args.mode+'\n') 
        
        if len(args.popfiles) > 0: 
            self.show('   Input Files: '+",".join([sf.name.split('/')[-1] for sf in args.popfiles])+'\n') 
        else: 
            self.show('   Input Files: '+",".join([sf.name.split('/')[-1] for sf in [args.prs,args.pheno]])+'\n') 
        self.show('Output Prefix: '+self.args.out+'\n\n')
        self.loc = None
        self.spl = '' 

    def initialize(self, f_name, t): 
        self.show('Beginning: '+f_name+'\n') 

    def show(self, msg, space=''):
        if space == 'NA': myspace = '' 
        else:             myspace = self.space
        if self.ACTIVE: 
            self.out2.write(myspace+msg)
            self.out2.flush() 

    def start(self,msg):
        if self.loc is not None: self.show('...Finished\n','NA') 
        self.space = '' 
        self.show('\n'+msg+':\n') 
        self.loc = None 
        self.space = '  '

    def end(self, msg,CHECK=[]): 
        self.spl += 3 
        if self.loc is not None: 
            self.show('...Finished ('+msg+')\n',space='NA') 
            if len(CHECK) == 2: 
                ws = ''.join([' ' for z in range(self.spl)]) 
                if CHECK[0] == 'VDATA': 
                    if abs(CHECK[1][0]) < 1000: self.show(ws+'Warning: Low sample size (<2k) may be unreliable\n') 
                    if abs(CHECK[1][1]) > 1 or abs(CHECK[1][1]) > 1: self.show(ws+'Warning: Skew/Kurtosis Suggests Non-Normality (try --normalize)\n') 
        self.loc = None 
        return

    def finish(self): 
        if self.loc is not None: self.show('...Finished\n','NA') 
        self.show('\nSibArc Completed\n','NA') 

    def update(self, msg): 
        if self.loc is not None: self.show('...Finished\n','NA') 
        self.loc = msg 
        self.show(msg+'...') 
        self.spl = len(msg)+3 























#############################################################################################
##################################  TRAIT PLOT     ##########################################
#############################################################################################



class PopPlot:
    def __init__(self, args, progress, results):
        self.args, self.progress, self.results = args, progress, results
        
        self.names = [k for k in self.results.keys()] 


        #if not self.args.normalize: self.fig_prefix = self.args.out+'.fig' 
        #else:                       self.fig_prefix = self.args.out+'.normalized.fig' 
        #if self.args.savePlotdata:
        #    self.out = open(self.args.out+'.plotData','w') 
        #    self.out.write('%-30s %40s %15s %s\n' % ('---', 'name', 'dataType', 'values')) 
        #self.names = args.names 
        #self.color_rainbow = plt.get_cmap('coolwarm')
        #self.herits = [round(i * 0.05, 2) for i in range(21)]
        #self.h_colors = [self.color_rainbow(1.*i/float(len(self.herits))) for i in range(len(self.herits))]
        #self.color_key = {h: c for h, c in zip(self.herits, self.h_colors)}
        #self.rounds = 1 
        #self.flen = min(6, len(self.names)) 
        #self.setup() 

    def make_curves(self): 
        self.fig = matplotlib.pyplot.gcf()
        self.axes, self.ax_index = [], 0  
        hl = int(1+len(self.names)/2.0) 
        if len(self.names) == 1: self.WD, self.HT, self.fs1, self.rows, self.cols = 10, 8, 15, 1, 1  
        elif len(self.names) == 2: self.WD, self.HT, self.fs1, self.rows, self.cols = 17, 8, 15, 1, 2  
        elif len(self.names) == 3: self.WD, self.HT, self.fs1, self.rows, self.cols = 24, 8, 15, 1, 3  
        elif len(self.names) == 4: self.WD, self.HT, self.fs1, self.rows, self.cols = 17, 15, 15, 2, 2  
        else:  self.WD, self.HT, self.fs1, self.rows, self.cols = 15, 15+(5*hl-2), 15, hl, 2  
        self.fig.set_size_inches(self.WD, self.HT) 
        for i in range(self.rows): 
            for j in range(self.cols):  
                self.axes.append(plt.subplot2grid((self.rows, self.cols), (i,j),     rowspan=1, colspan=1))
    
        for i,n in enumerate(self.names): 
            self.draw_curve(self.axes[i], n) 
        fig_name = self.args.out+'-dists'
        plt.subplots_adjust(left=0.1, bottom=0.15, right=0.95, top=0.94,wspace=0.2, hspace=0.2) 
        if  self.args.plotFormat == 'pdf': plt.savefig(fig_name+'.pdf', dpi=300) 
        else:                              plt.savefig(fig_name+'.png', dpi=300) 
        self.fig.clf() 

    def draw_curve(self, ax, n): 
        X, Y1, Y2 = self.results[n]['means'] 
        ax.plot(X, Y1, color = 'darkorange') 
        ax.scatter(X, Y2) 
        ax.set_title(n, fontsize=self.fs1) 
        ax.set_yticks([]) 
        ax.set_xlabel('Trait Percentile',fontsize=self.fs1) 
        ax.set_ylabel('PRS', fontsize=self.fs1) 
        #dict_keys(['params', 'pvals', 'effects', 'empirical', 'qc_val', 'QC', 'means'])


    def make_forests(self): 
        self.fig = matplotlib.pyplot.gcf()
        self.axes, self.ax_index = [], 0  
        hl = int(len(self.names)/4.0)
        self.WD, self.HT, self.fs1, self.rows, self.cols = 10, 3+(hl), 35, 1, 1  
        self.fig.set_size_inches(self.WD, self.HT) 
        self.ax = plt.subplot2grid((self.rows, self.cols), (0,0),     rowspan=1, colspan=1)
        cmap = plt.cm.get_cmap('hsv', len(self.names)) 


        self.yloc = 0 
        for i,n in enumerate(self.names): self.add_forest(n, cmap(i)) 
        self.yloc += 0.5
        self.ax.set_yticks([]) 
        
        self.ax.set_xticks([-2,-1.5,-1,-0.5,0.5,1,1.5,2]) 
        self.ax.set_xticklabels([-1*(-2+self.rO),-1*(-1.5+self.rO),-1*(-1+self.rO),-1*(-0.5+self.rO),0.5-self.rO,1-self.rO,1.5-self.rO,2-self.rO]) 

        self.ax.set_xlim(-2.2,2.2) 
        self.ax.set_ylim(self.yloc,1) 
        
        self.ax.set_xlabel('Tail Effects', fontsize=self.fs1) 


        for xl in [-2,-1.5,-1,-0.5,0.5,1,1.5,2]: 
            self.ax.plot([xl,xl],[self.yloc, 1], color='k', linestyle='--', alpha=0.1) 


        fig_name = self.args.out+'-effects'
        plt.subplots_adjust(left=0.05, bottom=0.25, right=0.95, top=0.94,wspace=0.2, hspace=0.2) 
        if  self.args.plotFormat == 'pdf': plt.savefig(fig_name+'.pdf', dpi=300) 
        else:                              plt.savefig(fig_name+'.png', dpi=300) 
        self.fig.clf() 



    def add_forest(self, n, clr): 

        self.rO =  1
        p1, p2 = self.results[n]['pvals'] 
        eL, eL1, eL2  = [-1*(self.rO + x) for x in self.results[n]['effects'][0]] 
        eH, eH1, eH2  = [self.rO + (x * -1) for x in self.results[n]['effects'][1]] 
        QC = self.results[n]['QC'] 
        self.ax.text(0, self.yloc, n, ha='center', va='center',fontsize=self.fs1) 
        self.ax.plot([eL1, eL2], [self.yloc, self.yloc], color = clr, lw = 2) 
        self.ax.plot([eH1, eH2], [self.yloc, self.yloc], color = clr, lw = 2) 
        self.ax.scatter(eL, self.yloc, color = clr, s = 200) 
        self.ax.scatter(eH, self.yloc, color = clr, s = 200)  
        self.yloc -= 1




        #dict_keys(['params', 'pvals', 'effects', 'empirical', 'qc_val', 'QC', 'means'])

    def setup(self):
        self.fig = matplotlib.pyplot.gcf()
        self.axes, self.ax_index = [], 0  
        xs1,xs2,xs3 = 16, 4,7
        ys1,ys2,ys3 = 16, 4,4

        self.fs0,self.fs1, self.fs2, self.fs3 = 15,14, 11.5, 9 
        

        cnt, row, col = 0, 0, 1 
        
        for i in range(self.rows): 
            for j in range(self.cols):  
                self.axes.append(plt.subplot2grid((self.rows, self.cols), (i,j),     rowspan=1, colspan=1))
        

        


    def add_popout(self): 
        ax = self.axes[self.ax_index] 

        b0, b1 = self.popOut.model.params

        X = self.popOut.phenoPrsMeans['pheno'] 
        Y = self.popOut.phenoPrsMeans['prs'] 
        Ye = [b0 + b1 * x for x in X] 
        ax.plot(self.popOut.qt_keys, Ye, color = 'darkorange') 
        ax.scatter(self.popOut.qt_keys, Y) 

        #args
        #progress
        #result
        #pair_data
        #ids
        #prs
        #phenos
        #prsPheno
        #phenoPrs
        #prsPhenoBins
            #phenoPrsBins
            #qt_keys
            #prsPhenoMeans
            #phenoPrsMeans
            #model
            #all_pvs
            #qc_pv
            #plot






    def draw_summary(self, f_name, name, pairs, h2, tt): 
        self.pairs, self.h2, self.tt = pairs, h2, tt 
        
        self.yMax, self.X = 0.5, sorted(self.pairs.keys())  
        self.m1 = [np.mean([p[0] for p in self.pairs[x]]) for x in self.X]
        self.m2 = [np.mean([p[1] for p in self.pairs[x]]) for x in self.X]
        self.s2 = [stats.sem([p[1] for p in self.pairs[x]]) for x in self.X]
        
        self.ax, ax2, ax3 = self.axes[self.ax_index], self.axes[self.ax_index + 1], self.axes[self.ax_index+2] 
        
        self.draw_herit_lines() 
        self.draw_curves(name, self.h2.body) 
        self.ax.text(-1, self.yMax*0.98, " ".join(name.split('_')), va='top', ha='left', fontsize=self.fs0+12) 
        if self.args.savePlotdata: 
            self.out.write('%-30s %40s %15s %s\n' % (f_name, name, 'X', ','.join([str(x) for x in self.X]))) 
            self.out.write('%-30s %40s %15s %s\n' % (f_name, name, 'm1', ','.join([str(round(x,4)) for x in self.m1]))) 
            self.out.write('%-30s %40s %15s %s\n' % (f_name, name, 'm2', ','.join([str(round(x,4)) for x in self.m2]))) 
            self.out.write('%-30s %40s %15s %s\n' % (f_name, name, 'sem', ','.join([str(round(x,4)) for x in self.s2]))) 
        
        t1 = SibTable(self.args, self.axes[self.ax_index+1]).make_res(tt,  self.mode, self.flen)  
        t2 = SibTable(self.args, self.axes[self.ax_index+2]).make_h2(h2,   self.mode, self.flen)  
        self.ax_index += 3


    def draw_herit_lines(self): 
        for i, h in enumerate(self.herits):
            sh = [(s * h)/2.0 for s in self.m1]
            if i == 0: self.ax.plot(self.X, sh, color=self.color_key[h], linewidth=3, alpha=0.7)
            else:
                z1 = [(a+b+b)/3.0 for a, b in zip(sh, lp)]
                z2 = [(a+b)/2.0 for a, b in zip(sh, lp)]
                self.ax.plot(self.X, z1, color=self.color_key[h], linewidth=4, alpha=0.3)
                self.ax.plot(self.X, z2, color=self.color_key[h], linewidth=4, alpha=0.3)
                self.ax.plot(self.X, sh, color=self.color_key[h], linewidth=4, alpha=0.3)
            lp = [z for z in sh]
            if max(lp) > self.yMax: self.yMax = max(lp) 





        
        


    def finish(self): 
        fig_name = self.args.out
        plt.subplots_adjust(left=0.1, bottom=0.15, right=0.95, top=0.94,wspace=0.2, hspace=0.2) 
        plt.savefig(fig_name+'.png', dpi=300) 
        plt.savefig(fig_name+'.pdf', dpi=300) 

    


    




class PopOut:
    def __init__(self, args, progress): 
        self.args, self.progress = args, progress 
        self.results = dd(lambda: {}) 
    
    def PopError(self, msg):  
        sys.stderr.write('POPError- '+msg+'\n') 
        sys.exit() 
    
    def set_pair_files(self,prsFile, phenoFile): 
        self.pair_data = dd(list) 
        for line in prsFile:
            line = line.split() 
            try: xId, xPrs = line[0], float(line[-1]) 
            except ValueError: continue 
            self.pair_data[xId] = [xPrs]  

        for line in phenoFile: 
            line = line.split() 
            try: xId, xVal = line[0], float(line[-1]) 
            except ValueError: continue 
            if xId in self.pair_data: 
                self.pair_data[xId].append(xVal) 

        self.ids, self.prs, self.phenos, self.prsPheno, self.phenoPrs = [], [], [], [], [] 
        for xId,xVals in self.pair_data.items(): 
            if len(xVals) == 2: 
                self.ids.append(xId) 
                self.prs.append(xVals[0]) 
                self.phenos.append(xVals[1]) 
                self.prsPheno.append([xVals[0], xVals[1]]) 
                self.phenoPrs.append([xVals[1], xVals[0]]) 
        self.summarize() 
        return self
   
    
    def process(self): 

        plot = PopPlot(self.args, self.progress, self.results) 
        plot.make_curves() 
        plot.make_forests() 
        self.write_result() 




    def write_result(self): 

        w = sys.stdout 
        w = open(self.args.out+'-result.txt','w') 
        w.write('%-18s %8s %9s %9s' % ('---', 'len', 'pv1','pv2')) 
        w.write(' %7s %7s %7s %7s %7s %7s' % ('e1', 'e1Lo', 'e1Hi','e2','e2Lo','e2Hi')) 
        w.write(' %6s %6s %6s' % ('QC', 'emp1','emp2')) 
        
        w.write('\n') 
        for r,K in self.results.items(): 
            rl = K['len'] 
            p1, p2 = K['pvals'] 
            e1, e2 = K['effects'] 
            emp, QC = K['empirical'], K['QC'] 
            w.write('%-18s %8d' % (r, K['len'])) 
            if p1 > 0.001: w.write(' %9.3f' % p1) 
            else:          w.write(' %9.1e' % p1) 
            if p2 > 0.001: w.write(' %9.3f' % p2) 
            else:          w.write(' %9.1e' % p2) 
            w.write(' %7.2f %7.2f %7.3f' % (e1[0], e1[1], e1[2])) 
            w.write(' %7.2f %7.2f %7.3f' % (e1[0], e1[1], e1[2])) 
            w.write(' %6s %6.2f %6.2f' % (QC, emp[0], emp[1])) 
            w.write('\n') 


    def read_and_run(self, f, name): 

        self.name, self.phenoPrs = name, [] 
        for i,lp in enumerate(f): 
            line = lp.split()
            try: pheno, prs = float(line[0]), float(line[1]) 
            except ValueError: 
                if i == 0: continue 
                else: self.PopError('Invalid Line: '+lp.strip()) 
            self.phenoPrs.append([pheno,prs]) 
        
        
        self.sampleLen = len(self.phenoPrs) 
        self.results[self.name]['len'] = len(self.phenoPrs) 

        self.phenoPrs.sort(key = lambda X: X[0]) 
        for i,(pheno,prs) in enumerate(self.phenoPrs): self.phenoPrs[i].append(round(100*(i/len(self.phenoPrs)),2))
        self.phenoPrs.sort(key = lambda X: X[1]) 
        for i,(pheno,prs, pheno_qt) in enumerate(self.phenoPrs): self.phenoPrs[i].append(round(100*(i/len(self.phenoPrs)),2))
        self.phenoPrs.sort() 
        self.phenos, self.prs = [d[0] for d in self.phenoPrs], [d[1] for d in self.phenoPrs] 
        
        self.run_test() 
        self.run_empirical_pv() 
        self.run_qc() 
        self.get_means() 



    def get_means(self): 
        prs_obs = dd(list) 
        prs_exp = dd(list) 
        yInt, beta = self.results[self.name]['params'] 
        for pheno, prs, pheno_qt, prs_qt in self.phenoPrs: 
            bin_val = int(pheno_qt) 
            prs_obs[bin_val].append(prs) 
            prs_exp[bin_val].append(yInt + beta*pheno) 
        X = [x for x in prs_obs.keys()] 
        prs_obs = [np.mean(prs_obs[x]) for x in X] 
        prs_exp = [np.mean(prs_exp[x]) for x in X] 
        self.results[self.name]['means'] = [X, prs_exp, prs_obs] 
        return 





    def assign_to_quants(self, y, n):
        q = np.percentile(y, np.linspace(100 / n, 100 * (n - 1) / n, n - 1))  # Quantiles at (1/n, 2/n, ..., (n-1)/n)
        y_cat = np.digitize(y, q) # Assign values to quantiles
        return y_cat





    def run_test(self, K=0.01):
        #prs = self.prs 
        #y = self.phenos
        # Fit the regression model
        X = sm.add_constant(self.phenos)  # Add a constant to the predictor variable (for the intercept)
        self.model = sm.OLS(self.prs, X).fit()

        self.results[self.name]['params'] = self.model.params 


        # Test upper tail "regression to mean"
        T_upper = np.percentile(self.phenos, 100 - K * 100)
        upper_tail_indices = np.where(self.phenos > T_upper)[0]
        residuals_upper = self.model.resid[upper_tail_indices]
        t_stat_upper, p_value_upper = stats.ttest_1samp(residuals_upper, 0)
        ci_upper = stats.t.interval(0.95, len(residuals_upper) - 1, loc=np.mean(residuals_upper), scale=np.std(residuals_upper) / np.sqrt(len(residuals_upper)))
        
        # Test lower tail "regression to mean"
        T_lower = np.percentile(self.phenos, K * 100)
        lower_tail_indices = np.where(self.phenos < T_lower)[0]
        residuals_lower = self.model.resid[lower_tail_indices]
        t_stat_lower, p_value_lower = stats.ttest_1samp(residuals_lower, 0)
        ci_lower = stats.t.interval(0.95, len(residuals_lower) - 1, loc=np.mean(residuals_lower), scale=np.std(residuals_lower) / np.sqrt(len(residuals_lower)))

        # Standard deviation of prs
        ms = np.std(self.prs)
        # Organizing the result
        
        self.results[self.name]['pvals'] = [p_value_lower, p_value_upper] 
        self.results[self.name]['effects'] = [[np.mean(residuals_lower) / ms, ci_lower[0] / ms, ci_lower[1]/ms]]
        self.results[self.name]['effects'].append([np.mean(residuals_upper) / ms, ci_upper[0] / ms, ci_upper[1]/ms]) 
        return


        my_result = {
            'effect.lower': np.mean(residuals_lower),
            'ci.lower.lower': ci_lower[0],
            'ci.lower.upper': ci_lower[1],
            'p.lower': p_value_lower,
            'effect.upper': np.mean(residuals_upper),
            'ci.upper.lower': ci_upper[0],
            'ci.upper.upper': ci_upper[1],
            'p.upper': p_value_upper,
            'sd': ms,
            're.lower': np.mean(residuals_lower) / ms,
            'ri.lower.lower': ci_lower[0] / ms,
            'ri.lower.upper': ci_lower[1] / ms,
            're.upper': np.mean(residuals_upper) / ms,
            'ri.upper.lower': ci_upper[0] / ms,
            'ri.upper.upper': ci_upper[1] / ms
        }
        for x,y in my_result.items(): self.result[x] = y 
        return



    def run_qc(self,stepsize=1): 
        trunc_data = [d for d in self.phenoPrs if d[2] > 10 and d[2] < 90] 
        X = sm.add_constant([t[0] for t in trunc_data])  # Add a constant to the predictor variable (for the intercept)
        self.model = sm.OLS([t[1] for t in trunc_data], X).fit()
        yInt, beta = self.model.params  
        resids = dd(list) 
        for pheno, prs, pheno_qt, prs_qt in trunc_data: 
            pred = yInt + beta*pheno 
            resids[int(pheno_qt)].append(pred-prs) 
        middle_pvs = [] 
        for k,R in resids.items(): 
            t_stat, pval = stats.ttest_1samp(R,0) 
            middle_pvs.append(pval) 

        qc_val = 80*min(middle_pvs) 
        self.results[self.name]['qc_val'] = qc_val 
        self.results[self.name]['QC'] = qc_val > 0.05 
        #self.qc_test_val = 80*min(middle_pvs) 





    def run_empirical_pv(self, n_q=100):
        q = self.assign_to_quants(self.phenos, n_q)
        self.all_pvs, stat_pairs = [], [] 
        # Loop through quantiles to perform t-tests
        for i in range(n_q):
            ptr = np.where(q == i)[0]
            residuals = self.model.resid[ptr]
            t_stat, p_value = stats.ttest_1samp(residuals, 0)
            stat_pairs.append([t_stat, i]) 
            self.all_pvs.append(p_value) 
        stat_pairs.sort() 
        stat_pairs = sorted([[sp[1],i] for i,sp in enumerate(sorted(stat_pairs))]) 
        p = [(100-stat_pairs[0][1])/100.0, (stat_pairs[-1][1]+1.0)/100.0]
        self.results[self.name]['empirical'] = p 
        self.empirical = p 
        #self.empirical = p 
        return 









def run_script(args, command_line, EXTEND=True):

    progress = PopProgress(args,command_line) 
    popOut = PopOut(args, progress) 

    if len(args.popfiles) > 0: 
        if len(args.names) != len(args.popfiles): 
            names = [p.name.split('/')[-1].split('.')[0] for p in args.popfiles] 
            if len(list(set(names))) == len(args.popfiles): args.names = names 
            else: args.names = ['trait'+str(i+1) for i,n in enumerate(args.popfiles)] 
        for i,f in enumerate(args.popfiles): popOut.read_and_run(f, args.names[i])
        popOut.process() 



    
    #popOut.run_qc_test() 
    else: 
        if args.prs and args.pheno: 
            print('yes') 
            popOut.set_pair_files(args.prs, args.pheno) 
    

            popOut.run_test() 
    
            popOut.plot() 



        for x,y in popOut.result.items(): 
            print(x,y) 


    sys.exit() 



if __name__ == '__main__':
    import argparse, sys
    usage = "usage: ./%prog [options] data_file"
    parser = argparse.ArgumentParser()
    parser.allow_abbrev=True
    
    parser.add_argument('popfiles',type=argparse.FileType('r'),nargs='+') 
    #parser.add_argument('pheno',type=argparse.FileType('r')) 
    #parser.add_argument('prs',type=argparse.FileType('r')) 
    parser.add_argument('--prs',type=argparse.FileType('r'), help='Phenotype File') 
    parser.add_argument('--pheno',type=argparse.FileType('r'), help='Phenotype File') 
    parser.add_argument("--names", nargs='+',default=[],type=str,help="Trait Name(s)")
    parser.add_argument("--out", type=str,default='out',help="Output Prefix")
    parser.add_argument("--plotFormat", type=str,default='pdf',help="Output Prefix")
    parser.add_argument("--mode", type=str,default='tailsOnly',help="Output Prefix")
    parser.add_argument("--normalize", action='store_true', default=False,help="Rank Inverse Normal Transform Data (For NonNormal Input)") 
    parser.add_argument("--silent", action='store_true', default=False,help="Suppress Output Stream") 
    args = parser.parse_args() 

    #run_script(args.prs, args.pheno, arg, ' '.join(sys.argv)) 
    run_script(args, ' '.join(sys.argv)) 






