import glob
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import pyrateless

from os import path

plt.rc('pgf',  texsystem='pdflatex')
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{lmodern}']
plt.rcParams["figure.figsize"] = (6,6)

def success_filter(df):
    return df.loc[df['success'] == 1.0, :]

def load(directory=None, multiplier=1, only_success=True, refresh=False):
    '''load df's from directory, average the columns and return a new df'''

    filename = directory + "-" + str(only_success) + ".csv"
    if not refresh:
        try:
            return pd.read_csv(filename)
        except FileNotFoundError:
            pass

    # aggregate all csv files into a single df
    pattern = path.join(directory, "*.csv")
    df = pd.concat(
        [pd.read_csv(f) for f in glob.glob(pattern)]
    )

    # optionally filter out decoding failures
    if only_success:
        df = success_filter(df)

    # if "status" in df:
    #     print()
    #     print(df['success'].sum())
    #     df['success'] = df.apply(lambda row: row['status'] == -4 or row['success'] == 1, axis=1)
    #     # df.loc[df['status'] == -4.0, :].loc[df['success'] ] = 1.0
    #     print(df['success'].sum())

    # average by overhead
    df = pd.DataFrame(
        [df[df['overhead']==x].mean() for x in df['overhead'].unique()]
    )

    # normalize
    df['reloverhead'] = df['overhead'] / df['num_inputs']
    df['decoding_additions'] /= df['num_inputs']
    df['decoding_multiplications'] /= df['num_inputs']
    df['rowadds'] /= df['num_inputs']
    df['rowmuls'] /= df['num_inputs']
    df['inactivations'] /= df['num_inputs']
    df['failure'] = 1-df[['success']]

    # add complexity
    df['complexity'] = df['decoding_additions']
    df['complexity'] += df['decoding_multiplications']
    df['complexity'] *= multiplier
    df['complexity'] += df['rowadds'] * 8
    df['complexity'] += df['rowmuls'] * 8
    df.sort_values(by=['reloverhead'], inplace=True)

    df.to_csv(filename)
    return df

def failure_plot(dfs=None, descriptors=None):
    plt.figure()
    ax = plt.subplot()
    popts = list()
    for df, descriptor in zip(dfs, descriptors):
        print("{}:".format(descriptor))
        x = df['reloverhead']
        y = df['failure']
        plt.semilogy(x, y, '.', label=descriptor)

        # fit a*exp(bx)
        df_opt = df
        df_opt = df_opt.loc[df_opt['failure'] < 1e-1, :]
        df_opt = df_opt.loc[df_opt['failure'] > 1e-3, :]

        x = df_opt['reloverhead']
        y = df_opt['failure']
        plt.semilogy(x, y, '.')

        x_opt = df_opt['reloverhead']
        y_opt = df_opt['failure']

        f = lambda t,a,b: a*np.exp(b*t)
        ff = lambda t,p: p[0]*np.exp(p[1]*t)

        p0 = None
        if popts:
            p0 = (popts[-1]["a"], popts[-1]["b"])
        try:
            popt, _ = sp.optimize.curve_fit(f, x_opt, y_opt, p0=p0, bounds=([1, -np.inf], [np.inf, -1]))
        except RuntimeError:
            continue

        a, b = popt
        popts.append({"a": a, "b": b, "num_inputs": df['num_inputs'].mean()})
        t = np.linspace(0, x_opt.max())
        # plt.semilogy(t, ff(t,popt))

        f = lambda t,a,b: a*np.exp(b*t)
        n = df["num_inputs"].mean()
        # am = np.exp(0.005161818711756707*n)
        bm = -594.2043276421593
        am = 10*np.exp(4.06755643e-03*n)
        # plt.semilogy(t, f(t, am, bm))

        # plt.semilogy(x, y-f(x, am, bm), '.', label=descriptor)



    # plot binary random fountain
    overhead = np.arange(0, 10)
    plt.semilogy(overhead/df['num_inputs'].mean(), 1/np.power(2, overhead), label="$2^{-\delta}$")

    # plot lower bound
    # K = 4000
    # M = 3998
    # delta = 0.9999999701976676
    # soliton = pyrateless.Soliton(delta=delta, symbols=K, mode=M)
    # overheads = np.linspace(0.25, 0.4, 10)
    # failure = [
    #     pyrateless.optimize.decoding_failure_prob_estimate(
    #         soliton=soliton,
    #         num_inputs=K,
    #         overhead=1+x,
    #     ) for x in overheads
    # ]
    # plt.semilogy(overheads, failure, label="bound")

    plt.setp(ax.get_xticklabels(), fontsize=25)
    plt.setp(ax.get_yticklabels(), fontsize=25)
    plt.xlabel("Overhead", fontsize=25)
    plt.ylabel("Failure Probability", fontsize=25)
    plt.ylim(ymax=1, ymin=1e-6)
    plt.xlim(xmin=0, xmax=0.1)
    # plt.xlim([0, 0.5])
    plt.grid()
    plt.legend()
    # plt.autoscale(enable=True)
    # plt.tight_layout()
    # plt.savefig("failure.png")

    # plot optimized parameters
    # filter out noise
    popts = pd.DataFrame(popts)
    popts = popts.loc[popts["a"] > 0, :]
    # popts = popts.loc[popts["a"] < 1e10, :]
    popts = popts.loc[popts["b"] < 0, :]
    # popts = popts.loc[popts["num_inputs"] <= 4000, :]

    # fit curve to parameter a
    x_opt = popts["num_inputs"]
    y_opt = popts["a"]
    fa = lambda t,a,b: a*np.exp(b*t)
    ffa = lambda t,p: p[0]*np.exp(p[1]*t)
    a0 = 1
    b0 = np.log(y_opt)
    b0 /= x_opt
    b0 = b0.mean()
    print("(a0, b0)=({}, {})".format(a0, b0))
    popt, _ = sp.optimize.curve_fit(
        fa, x_opt, y_opt, p0=(a0, b0), bounds=(0, [1e2, np.inf]),
    )

    print("meta-fit parameters to a: {}".format(popt))

    plt.figure()
    ax = plt.subplot()
    t = np.linspace(x_opt.min(), x_opt.max())
    plt.semilogy(t, ffa(t, popt))
    plt.semilogy(t, fa(t, a0, b0), '--')
    plt.semilogy(x_opt, y_opt, '.')
    plt.title("$a$")

    plt.figure()
    ax = plt.subplot()
    x_opt = popts["num_inputs"]
    y_opt = popts["b"]
    c0 = y_opt.mean()
    print("meta-fit parameters to b: {}".format(c0))
    t = np.linspace(x_opt.min(), x_opt.max())
    plt.plot(t, c0*np.ones(len(t)))
    plt.plot(popts["num_inputs"], popts["b"], '.')
    plt.title("$b$")

    return

def required_overhead_plot(dfs=None, descriptors=None, tfps=[1e-1, 1e-2, 1e-3]):
    r = list()
    for df in dfs:
        dct = dict()
        dct["x"] = df["num_inputs"].mean()
        for tfp in tfps:
            i = (df["failure"]-tfp).abs().idxmin(axis=1)
            dct[str(tfp)+"-actual"] = df["failure"][i]
            dct[str(tfp)] = df["reloverhead"][i]
            dct[str(tfp)+"-complexity"] = df["complexity"][i]

        r.append(dct)

    df = pd.DataFrame(r)
    print(df.head())
    plt.figure()
    ax = plt.subplot()
    for tfp in tfps:
        plt.plot(df["x"], df[str(tfp)], '.', label="tfp={}".format(tfp))

    plt.grid()
    plt.title("required overhead")
    plt.legend()

    plt.figure()
    ax = plt.subplot()
    for tfp in tfps:
        plt.semilogy(df["x"], np.abs(tfp-df[str(tfp)+"-actual"]), '.', label="tfp={}".format(tfp))

    plt.grid()
    plt.title("tfp error")
    plt.legend()

    plt.figure()
    ax = plt.subplot()
    for tfp in tfps:
        plt.plot(df["x"], df[str(tfp)+"-complexity"], '.', label="tfp={}".format(tfp))

    plt.grid()
    plt.title("complexity")
    plt.legend()

    plt.figure()
    ax = plt.subplot()
    plt.plot(df["x"], (df[str(1e-2)] / df[str(1e-1)] - 1), '.', label="1e-1/1e-2")
    plt.plot(df["x"], (df[str(1e-3)] / df[str(1e-2)] - 1), '.', label="1e-2/1e-3")
    plt.plot(df["x"], (df[str(1e-3)] / df[str(1e-1)] - 1)/2, '.', label="1e-1/1e-3")
    plt.grid()
    plt.title("relative")
    plt.legend()

    plt.figure()
    ax = plt.subplot()
    plt.plot(df["x"], (df[str(1e-2)] - df[str(1e-1)])/1, '.', label="1e-1, 1e-2")
    plt.plot(df["x"], (df[str(1e-3)] - df[str(1e-2)])/1, '.', label="1e-2, 1e-3")
    plt.plot(df["x"], (df[str(1e-3)] - df[str(1e-1)])/2, '.', label="1e-2, 1e-3")
    plt.grid()
    plt.title("absolute")
    plt.legend()

    # plt.xlabel("num_inputs")
    # plt.ylabel("(required relative overhead)")
    return

def complexity_plot(dfs=None, descriptors=None, complexity_coefficients=None):
    if complexity_coefficients is None:
        complexity_coefficients = np.ones(len(dfs))
    plt.figure()
    ax = plt.subplot()
    for df, descriptor, multiplier in zip(dfs, descriptors, complexity_multipliers):
        print("{}:".format(descriptor))
        # print(df.head())
        plt.semilogy(df['reloverhead'], df['complexity'], label=descriptor)

    plt.setp(ax.get_xticklabels(), fontsize=25)
    plt.setp(ax.get_yticklabels(), fontsize=25)
    plt.xlabel("Overhead", fontsize=25)
    plt.ylabel("Complexity", fontsize=25)
    # plt.xlim([0, 0.5])
    plt.grid()
    plt.legend()
    plt.autoscale(enable=True)
    # plt.tight_layout()
    # plt.savefig("complexity.png")
    return

# R10 failure-overhead curve approximation
# R10(1000), 1e-1, 1e-2 at 1.5%, 2%
# R10(2000), 1e-1, 1e-2 at 2.1%, 2.4%
# R10(3000), 1e-1, 1e-2 at 2.9%, 3.3%
# R10(4000), 1e-1, 1e-2 at 3.7%, 4.1%
# R10(6000), 1e-1, 1e-2 at 5.4%, 5.6%
# R10(4000), 1e-1, 1e-2 at 7%, 7.5%


if __name__ == '__main__':
    refresh = False
    prefix = "./simulations/"
    # directories = [d for d in glob.glob("./simulations/R10*") if path.isdir(d)]
    directories = [
        "R10(1000)",
        "R10(1100)",
        "R10(1200)",
        "R10(1300)",
        "R10(1400)",
        "R10(1500)",
        "R10(1600)",
        "R10(1700)",
        "R10(1800)",
        "R10(1900)",
        "R10(2000)",
        "R10(2100)",
        "R10(2200)",
        "R10(2300)",
        "R10(2400)",
        "R10(2500)",
        "R10(2600)",
        "R10(2700)",
        "R10(2800)",
        "R10(2900)",
        "R10(3000)",
        "R10(3100)",
        "R10(3200)",
        "R10(3300)",
        "R10(3400)",
        "R10(3500)",
        "R10(3600)",
        "R10(3700)",
        "R10(3800)",
        "R10(4000)",
        "R10(6000)",
        "R10(7000)",
    ]
    directories = [
        "R10(1000)",
        "R10(2000)",
        "R10(3000)",
        "R10(4000)",
        "R10(6000)",
        "R10(7000)",
        "R10(8000)",
    ]
    directories = [prefix + d for d in directories]
    descriptors = [d.strip(prefix) for d in directories]
    print(directories)
    print(descriptors)

    # directories = [
    #     "R10Parameters(1000)",
    #     "R10(2000)",
    #     "R10(3000)",
    #     "R10Parameters(4000)",
    #     "R10(6000)",
    #     "R10(8000)",
    #     # "LTParameters(4000, Soliton(4000, 3998, 0.9999999701976676))",
    #     # "QLTParameters{DT,UInt8}(4000, Soliton(4000, 3998, 0.9999999701976676))",
    #     # "LTParameters(1000, Soliton(1000, 998, 0.9999999701976676))",
    #     # "QLTParameters{DT,UInt8}(1000, Soliton(1000, 998, 0.9999999701976676))",
    #     # "LTQ{Float64,DT}(4000, Soliton(4000, 3998, 0.9999999701976676))",
    #     # "LTQ{Float64,DT}(4000, Soliton(4000, 3998, 0.9999999701976676))-1",
    #     # "LTParameters(4000, Soliton(4000, 3998, 0.999))",
    #     # "QLTParameters{DT,UInt8}(4000, Soliton(4000, 3998, 0.999))",
    # ]
    # directories = [prefix + directory for directory in directories]
    # descriptors = [
    #     "R10(1000)",
    #     "R10(2000)",
    #     "R10(3000)",
    #     "R10(4000)",
    #     "R10(6000)",
    #     "R10(8000)",
    #     # "LT(4000, 3998, 0.9999...)",
    #     # "LT256(4000, 3998, 0.9999...)",
    #     # "LT(1000, 998, 0.9999...)",
    #     # "LT256(1000, 998, 0.9999...)",
    #     # "LT-Float64(4000, 3998, 0.9999...)",
    #     # "LT-Float64-1e10(4000, 3998, 0.9999...)",
    # ]
    # complexity_multipliers = [1, 1, 1, 1, 1, 1, 1, 8, 1, 8, 8*4, 8*4]
    complexity_multipliers = np.ones(len(directories))

    # aggregate dataframes
    dfs_success = [load(f, only_success=True, refresh=refresh) for f in directories]
    dfs_any = [load(f, only_success=False, refresh=refresh) for f in directories]

    # make plots
    failure_plot(dfs_any, descriptors)
    # complexity_plot(dfs_success, descriptors, complexity_multipliers)
    # required_overhead_plot(dfs_any, descriptors)
    plt.show()
