"""
Independent re-implementation of the ER-027 Appendix A quartering simulation.

Written separately from cpa_quarter_sim.py, from the parameters stated in the
appendix only (absorption half-life 0.90 h, distribution-phase half-life 2.9 h
carrying 34% of disposition AUC, terminal half-life 38.21 h; symmetric Dirichlet
quarter shares; tablet content CV 2%; 120 days, days 1-30 discarded), to check
that the published results do not depend on one implementation.

It uses 2,000 simulated patients per condition and its own random seed, so it
matches the published results to about one percentage point, not exactly. See
REPRODUCTION.md.

Run:  python check_independent.py
"""
import numpy as np
ln2=np.log(2)
ka=ln2/0.90; al=ln2/2.9; be=ln2/38.21
# choose A,B so that A/al : B/be = 0.34 : 0.66 (distribution phase share of AUC from the disposition terms)
B=1.0; A=(0.34/0.66)*B*al/be
def C(t):  # single unit-dose curve, t in h (array), zero before 0
    t=np.asarray(t); out=A*np.exp(-al*t)+B*np.exp(-be*t)-(A+B)*np.exp(-ka*t); out[t<0]=0; return out
# check tmax and auc share
tt=np.arange(0,400,0.05); c=C(tt); print('tmax',tt[c.argmax()].round(2),'AUC fast share',round((A/al-0)/(A/al+B/be),2))
days=120; K=np.arange(days*24*2)/2.0
kernel=C(np.arange(0,days*24,0.5))
rng=np.random.default_rng(1)
def troughs(doses):
    # doses per day (relative), compute conc on grid via convolution
    dgrid=np.zeros(len(K)); dgrid[np.arange(days)*48]=doses
    conc=np.convolve(dgrid,kernel)[:len(K)]
    tr=conc[np.arange(31,days)*48-1]  # just before dose
    return tr
even=troughs(np.ones(days))
def dirichlet_cv(cv,n):
    # symmetric Dirichlet with 4 parts: CV of share = sqrt((1-p)/(p*(a0+1))) with p=1/4
    a0=(3/(cv**2))-1; return rng.dirichlet([a0/4]*4,size=n)*4
def run(scen,cv,N=2000,miss=False,eod=False):
    out=[]
    for i in range(N):
        if scen=='home':
            q=dirichlet_cv(cv,days//4)
            for row in q: rng.shuffle(row)
            d=q.ravel()[:days]*rng.normal(1,0.02,days//4).repeat(4)[:days]
        elif scen=='pool':
            d=dirichlet_cv(cv,days)[:,0]*rng.normal(1,0.02,days)
        elif scen=='miss':
            d=np.ones(days)
            for s in range(0,days,28):
                L=min(28,days-s)
                if L==28 or rng.random()<L/28:  # partial last block: miss with probability L/28 (v1.1)
                    d[s+rng.integers(0,L)]=0
        elif scen=='eod':
            d=np.tile([2,0],days//2).astype(float)
        out.append(troughs(d)/even)
    o=np.concatenate(out); return np.percentile(o,[2.5,97.5]).round(3), (o<0.8).mean().round(4)
for s,cv in [('home',.15),('pool',.15),('home',.30),('pool',.30)]: print(s,cv,run(s,cv))
print('miss',run('miss',0)); print('eod',run('eod',0,N=10))
print('accum',1/(1-np.exp(-24*be)))
