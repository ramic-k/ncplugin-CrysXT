"""Does isotropic sampling recover the powder average? (the natural first question)

Two checks, both should recover the random powder:
  1. R -> 1 limit (texture off): the textured cross section equals the random-powder xs.
  2. Isotropic incidence: averaging the ORIENTED sampler over isotropically-distributed incident
     directions reproduces the powder scattering-angle distribution -- i.e. the texture conserves
     the angle-integrated scattering, at the sampling level (not just the cross section).

Run from the plugin root with the plugin loaded via NCRYSTAL_PLUGIN_LIST.  Material: Fe_sg229_CrysXT.
"""
import numpy as np, tempfile
import NCrystal as NC

NCMAT = "data/Fe_sg229_CrysXT.ncmat"

def strip_custom(path):
    out, inc = [], False
    for ln in open(path):
        s = ln.strip()
        if s.startswith("@CUSTOM_CRYSXT"): inc = True; continue
        if inc and s.startswith("@"): inc = False
        if inc: continue
        out.append(ln)
    return "".join(out)

base = strip_custom(NCMAT)
def w(custom):
    f = tempfile.NamedTemporaryFile("w", suffix=".ncmat", delete=False)
    f.write(base + ("\n" + custom if custom else "")); f.close(); return f.name

# ---- 1) R -> 1 limit recovers the random powder cross section ----
TEX_R1 = "@CUSTOM_CRYSXT\n  Texture 0 0 1 1.0 0.5\n  Texture 0 0 1 1.0 0.5"
powder = NC.load(f"{w('')};dcutoff=0.3;comp=coh_elas")
texR1  = NC.load(f"{w(TEX_R1)};dcutoff=0.3;comp=coh_elas")
print("[1] R->1 textured xs vs random powder xs:")
ok1 = True
for wl in [1.5, 2.5, 3.5]:
    a, b = powder.scatter.xsect(wl=wl), texR1.scatter.xsect(wl=wl)
    rel = abs(b-a)/a*100; ok1 &= rel < 1e-3
    print(f"    wl={wl} A : powder={a:.4f}  tex(R=1)={b:.4f}  rel={rel:.4f}%")

# ---- 2) oriented sampler averaged over ISOTROPIC incidence -> powder scattering-angle dist ----
TEX_ORI = "@CUSTOM_CRYSXT\n  Texture 0.6 0.1 0.78 0.30 0.6\n  Texture 0.2 0.4 0.9 0.25 0.4"
ori = NC.load(f"{w(TEX_ORI)};dcutoff=0.3;comp=coh_elas;mos=0.0005deg;"
              "dir1=@crys_hkl:1,0,0@lab:1,0,0;dir2=@crys_hkl:0,1,0@lab:0,1,0")
wl = 2.8; ek = NC.wl2ekin(wl); rng = np.random.default_rng(0)
ep, dp = powder.scatter.sampleScatter(ek, (0,0,1), repeat=80000)
dp = np.asarray(dp); dp = dp.T if dp.shape[0]==3 else dp
mu_pw = dp @ np.array([0,0,1.])
ins = rng.normal(size=(1600,3)); ins /= np.linalg.norm(ins,axis=1)[:,None]; mu_or = []
for k in range(1600):
    e2, d2 = ori.scatter.sampleScatter(ek, tuple(ins[k]), repeat=50)
    d2 = np.asarray(d2); d2 = d2.T if d2.shape[0]==3 else d2
    mu_or.append(d2 @ ins[k])
mu_or = np.concatenate(mu_or)
bins = np.linspace(-1,1,25)
hp,_ = np.histogram(mu_pw, bins, density=True); ho,_ = np.histogram(mu_or, bins, density=True)
corr = np.corrcoef(hp,ho)[0,1]; maxd = np.max(np.abs(hp-ho))/np.max(hp)*100
ok2 = corr > 0.999 and maxd < 2.0
print("[2] isotropic-incidence oriented sampler vs powder (scattering-cosine distribution):")
print(f"    histogram correlation={corr:.4f}  max bin diff={maxd:.1f}%  mean cos: powder={mu_pw.mean():.4f} oriented-iso-avg={mu_or.mean():.4f}")
print(f"\nPOWDER RECOVERY: {'PASS' if (ok1 and ok2) else 'FAIL'}")
