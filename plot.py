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
plt.rcParams["figure.dpi"] = 200

def success_filter(df):
    return df.loc[df['success'] == 1.0, :]

def dfgen(pattern, index, only_success=False):
    for df in (pd.read_csv(f) for f in glob.glob(pattern)):
        if only_success:
            df = success_filter(df)
        df['count'] = 1
        df.set_index(index)
        df = df.groupby(index).sum()
        yield df
    return

def load(directory=None, index='overhead', only_success=True, refresh=False):
    '''load df's from directory, average the columns and return a new df'''
    filename = directory + "-" + str(only_success) + ".csv"
    print('loading', filename)
    if not refresh:
        try:
            return pd.read_csv(filename)
        except FileNotFoundError:
            pass
    pattern = path.join(directory, "*.csv")
    g = dfgen(pattern, index) # load dfs iteratively
    agg = next(g)
    for df in g:
        agg = agg.add(df, fill_value=0)

    # compute average
    for column in df:
        if column == 'overhead':
            continue
        df[column] /= df['count']

    # give the x variable a predictable name
    df['x'] = df[index]
    df.sort_values(by=['x'], inplace=True)

    # normalize by number of inputs
    for column in df:
        if column == 'num_inputs':
            continue
        df[column] /= df['num_inputs']

    # compute complexity
    df['complexity'] = df['decoding_additions']
    df['complexity'] += df['decoding_multiplications']
    df['complexity'] *= multiplier
    df['complexity'] += df['rowadds'] * 8
    df['complexity'] += df['rowmuls'] * 8

    df.to_csv(filename)
    return df

def load_old(directory=None, multiplier=1, only_success=True, refresh=False):
    '''load df's from directory, average the columns and return a new df'''
    filename = directory + "-" + str(only_success) + ".csv"
    print('loading', filename)
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

    # average by overhead or erasures
    if 'overhead' in df:
        count = [
            len(df[(df['overhead']==x)]) for x in df['overhead'].unique()
        ]
        failure_count = [
            len(df[(df['overhead']==x) & (df['success']==0)]) for x in df['overhead'].unique()
        ]
        df = pd.DataFrame(
            [df[df['overhead']==x].mean() for x in df['overhead'].unique()]
        )
        df['x'] = df['overhead'] / df['num_inputs']
        df['count'] = count
        df['failure_count'] = failure_count
    elif 'erasures' in df:
        count = [
            len(df[(df['erasures']==x)]) for x in df['erasures'].unique()
        ]
        failure_count = [
            len(df[(df['erasures']==x) & (df['success']==0)]) for x in df['erasures'].unique()
        ]
        print('fails', failure_count)
        df = pd.DataFrame(
            [df[df['erasures']==x].mean() for x in df['erasures'].unique()]
        )
        df['x'] = df['erasures'] / df['length']
        df['count'] = count
        df['failure_count'] = failure_count
    else:
        raise ValueError('df must have column overhead or erasures')

    # normalize
    df['decoding_additions'] /= df['num_inputs']
    df['decoding_multiplications'] /= df['num_inputs']
    df['rowadds'] /= df['num_inputs']
    df['rowmuls'] /= df['num_inputs']
    df['inactivations'] /= df['num_inputs']
    df['failure'] = 1-df[['success']]

    # compute complexity
    df['complexity'] = df['decoding_additions']
    df['complexity'] += df['decoding_multiplications']
    df['complexity'] *= multiplier
    df['complexity'] += df['rowadds'] * 8
    df['complexity'] += df['rowmuls'] * 8
    df.sort_values(by=['x'], inplace=True)

    df.to_csv(filename)
    return df

def failure_plot(dfs=None, descriptors=None):
    plt.figure()
    ax = plt.subplot()
    for df, descriptor in zip(dfs, descriptors):
        print("{}:".format(descriptor))
        print(df)
        df = df.loc[df['failure'] > 0, :] # remove points with no failures
        x = df['x']
        y = df['failure']
        plt.semilogy(x, y, '-o', label=descriptor)

    # plot binary random fountain
    # overhead = np.arange(0, 10)
    # plt.semilogy(overhead/df['num_inputs'].mean(), 1/np.power(2, overhead), label="$2^{-\delta}$")

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
    # plt.semilogy(overheads, failure, label="Bound, 4000")

    # K = 8000
    # M = K-2
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
    # plt.semilogy(overheads, failure, label="Bound, 8000")

    plt.setp(ax.get_xticklabels())
    plt.setp(ax.get_yticklabels())
    plt.xlabel("Erasure Rate")
    # plt.xlabel("Relative Overhead")
    plt.ylabel("FER")
    # plt.xlim((0.4, 0.5))
    # plt.ylim((1e-6, 1))
    # plt.xlim((0, 0.1))
    # plt.ylim((1e-4, 1))
    plt.legend()
    plt.title("LDPC (2400, 4800) Failure-Erasure Curve")
    # plt.title("Raptor Failure-Overhead Curve")
    # plt.autoscale(enable=True)
    plt.tight_layout()
    # plt.savefig("ldpc-failure-overhead.png", dpi="figure")
    # plt.savefig("./plots/180419/rq.png", dpi="figure")
    return

def required_overhead_plot(dfs=None, descriptors=None, tfps=[1e-1, 1e-2, 1e-3]):
    r = list()
    for df in dfs:
        dct = dict()
        dct["x"] = df["num_inputs"].mean()
        for tfp in tfps:
            i = (df["failure"]-tfp).abs().idxmin(axis=1)
            dct[str(tfp)+"-actual"] = df["failure"][i]
            dct[str(tfp)] = df["x"][i]
            dct[str(tfp)+"-complexity"] = df["complexity"][i]

        r.append(dct)

    df = pd.DataFrame(r)
    print("Failure-Overhead")
    print(df.head())
    plt.figure()
    ax = plt.subplot()
    for tfp in tfps:
        x, y = df["x"], df[str(tfp)]
        plt.plot(x, y, '.', label="$P_F={}$".format(tfp))
        f = lambda t,a,b: a+b*t
        ff = lambda t,p: p[0]+p[1]*t
        popt, _ = sp.optimize.curve_fit(f, x, y)
        print("found parameters {} for fp={}".format(popt, tfp))
        t = np.linspace(x.min(), x.max())
        plt.plot(t, ff(t, popt))

    plt.xlim((1000, 8000))
    plt.ylim((0.01, 0.08))
    plt.title("Overhead at Given Failure Probability")
    plt.xlabel("\# Source Symbols")
    plt.ylabel("Relative Overhead")
    plt.legend()
    plt.tight_layout()
    plt.savefig("r10_required_overhead.png", dpi="figure")
    return

    plt.figure()
    ax = plt.subplot()
    for tfp in tfps:
        plt.semilogy(df["x"], np.abs(tfp-df[str(tfp)+"-actual"])/tfp, '.', label="tfp={}".format(tfp))

    plt.title("tfp error")
    plt.legend()

    print("Complexity")
    plt.figure()
    ax = plt.subplot()
    for tfp in tfps:
        x, y = df["x"], df[str(tfp)+"-complexity"]
        plt.plot(x, y, '.', label="tfp={}".format(tfp))
        f = lambda t,a,b: a+b*t
        ff = lambda t,p: p[0]+p[1]*t
        popt, _ = sp.optimize.curve_fit(f, x, y)
        print("found parameters {} for fp={}".format(popt, tfp))
        t = np.linspace(x.min(), x.max())
        plt.plot(t, ff(t, popt))

    plt.xlim((1000, 8000))
    plt.ylim((150, 550))
    plt.title("Complexity at Given Failure Probability")
    plt.xlabel("\# Source Symbols")
    plt.ylabel("Decoding Complexity")
    plt.legend()
    plt.tight_layout()
    # plt.savefig("r10_complexity.png", dpi="figure")
    return

    plt.figure()
    ax = plt.subplot()
    plt.plot(df["x"], (df[str(1e-2)] / df[str(1e-1)] - 1), '.', label="1e-1/1e-2")
    plt.plot(df["x"], (df[str(1e-3)] / df[str(1e-2)] - 1), '.', label="1e-2/1e-3")
    plt.plot(df["x"], (df[str(1e-3)] / df[str(1e-1)] - 1)/2, '.', label="1e-1/1e-3")
    plt.title("relative")
    plt.legend()

    plt.figure()
    ax = plt.subplot()
    plt.plot(df["x"], (df[str(1e-2)] - df[str(1e-1)])/1, '.', label="1e-1, 1e-2")
    plt.plot(df["x"], (df[str(1e-3)] - df[str(1e-2)])/1, '.', label="1e-2, 1e-3")
    plt.plot(df["x"], (df[str(1e-3)] - df[str(1e-1)])/2, '.', label="1e-2, 1e-3")
    plt.title("absolute")
    plt.legend()
    return

def required_complexity_plot(dfs, descriptors):
    pf = [1e-1, 1e-2, 1e-3]
    a = [5.43476844e-03, 9.24066372e-03, 1.36168053e-02]

    # average of 1e-1 and 1e-2 values since the distance to the measured pf is
    # much smaller for these points.
    b = (8.09230267e-06 + 8.01977332e-06) / 2
    f = lambda t,a,b: a+b*t

    r = list()
    for df, descriptor in zip(dfs, descriptors):
        K = df["num_inputs"].mean()
        x1 = f(K, a[0], b) # 1e-1
        x2 = f(K, a[1], b) # 1e-2
        x3 = f(K, a[2], b) # 1e-3

        dct = dict()
        dct["x"] = df['num_inputs'].mean()

        i = (df["x"]-x1).abs().idxmin(axis=1)
        dct["1e-1"] = df["complexity"][i]

        i = (df["x"]-x2).abs().idxmin(axis=1)
        dct["1e-2"] = df["complexity"][i]

        i = (df["x"]-x3).abs().idxmin(axis=1)
        dct["1e-3"] = df["complexity"][i]
        r.append(dct)

    df = pd.DataFrame(r)
    print(df)
    plt.figure()
    ax = plt.subplot()
    for tfp in ["1e-1", "1e-2", "1e-3"]:
        x, y = df["x"], df[tfp]
        plt.plot(x, y, '.', label="$P_F={}$".format(tfp))
        f = lambda t,a,b: a+b*t
        ff = lambda t,p: p[0]+p[1]*t
        popt, _ = sp.optimize.curve_fit(f, x, y)
        print("found parameters {} for fp={}".format(popt, tfp))
        t = np.linspace(x.min(), x.max())
        plt.plot(t, ff(t, popt))

        # x = [0, x1-eps, x1, x2, x3, x4]
        # y = [0, 0, 1-1e-1, 1-1e-2, 1-1e-3, 1]
        # l = ax.plot(x, y, label=descriptor)
        # color = l[-1].get_color()
        # x = df["x"]
        # y = np.array(df["complexity"])
        # plt.plot(x, y, '--', color=color)

    plt.title("R10 Complexity")
    plt.xlabel("\# Source Symbols")
    plt.ylabel("Complexity")
    plt.xlim((1000, 8000))
    plt.ylim((2e2, 6e2))
    plt.legend()
    plt.tight_layout()
    plt.savefig("r10_required_complexity.png", dpi='figure')
    return

def cdf_plot(dfs, descriptors):
    pf = [1e-1, 1e-2, 1e-3]
    a = [5.43476844e-03, 9.24066372e-03, 1.36168053e-02]

    # average of 1e-1 and 1e-2 values since the distance to the measured pf is
    # much smaller for these points.
    b = (8.09230267e-06 + 8.01977332e-06) / 2
    f = lambda t,a,b: a+b*t
    eps = np.finfo(float).eps
    plt.figure()
    ax = plt.subplot()
    for df, descriptor in zip(dfs, descriptors):
        K = df["num_inputs"].mean()
        x1 = f(K, a[0], b)
        x2 = f(K, a[1], b)
        x3 = f(K, a[2], b)
        x = [0, x1-eps, x1, x2, x3, 2*x3]
        y = [0, 0, 1-1e-1, 1-1e-2, 1-1e-3, 1]
        l = ax.plot(x, y, label=descriptor)
        color = l[-1].get_color()
        x = df["x"]
        y = np.array(df["success"])
        plt.plot(x, y, '--', color=color)

    plt.title("R10 Completion CDF")
    plt.xlabel("Relative Overhead")
    plt.ylabel("Probability")
    plt.xlim((0, 0.1))
    plt.ylim((0, 1))
    plt.legend()
    plt.tight_layout()
    plt.savefig("cdf.png", dpi='figure')
    return

def complexity_prediction_plot(dfs, descriptors):
    pf = [1e-1, 1e-2, 1e-3]
    a = [2.22413913e+02, 1.95991293e+02, 1.68044631e+02]
    b = [3.64861777e-02, 2.97971578e-02, 2.70341188e-02]
    f = lambda t,a,b: a+b*t
    plt.figure()
    ax = plt.subplot()
    for df, descriptor in zip(dfs, descriptors):
        K = df["num_inputs"].mean()
        y0 = 0
        y1 = f(K, a[0], b[0])
        y2 = f(K, a[1], b[1])
        y3 = f(K, a[2], b[2])
        x0 = 0
        x1 = 1e-1
        x2 = x1+1e-2
        x3 = x2+1e-3
        x = [x0, x1, x2, x3]
        y = [y0, y1, y2, y3]
        l = ax.plot(x, y, label=descriptor)
        color = l[-1].get_color()
        x = df["success"]
        y = np.array(df["complexity"])
        plt.plot(x, y, '--', color=color)

    plt.title("R10 Complexity")
    plt.xlabel("Failure Probability")
    plt.ylabel("Complexity")
    plt.xlim((0, 0.1))
    plt.ylim((0, 1))
    plt.legend()
    plt.tight_layout()
    plt.savefig("cdf.png", dpi='figure')
    return

def complexity_plot(dfs=None, descriptors=None, complexity_coefficients=None):
    if complexity_coefficients is None:
        complexity_coefficients = np.ones(len(dfs))
    plt.figure()
    ax = plt.subplot()
    for df, descriptor, multiplier in zip(dfs, descriptors, complexity_coefficients):
        print("{}:".format(descriptor))
        # print(df.head())
        plt.semilogy(df['x'], df['complexity'], label=descriptor)

    plt.title("Raptor Complexity")
    plt.xlabel("Relative Overhead")
    plt.ylabel("Complexity")
    plt.xlim((0, 0.1))
    plt.ylim((10, 2e3))
    plt.legend()
    # plt.autoscale(enable=True)
    # plt.tight_layout()
    plt.savefig("./plots/180419/raptor_complexity.png")
    return


def check_status(directory=None, multiplier=1, only_success=True, refresh=True):
    '''load df's from directory, average the columns and return a new df'''

    filename = directory + "-" + str(only_success) + ".csv"
    if not refresh:
        try:
            return pd.read_csv(filename)
        except FileNotFoundError:
            pass

    pattern = path.join(directory, "*.csv")
    for f in glob.glob(pattern):
        df = pd.read_csv(f)
        if "status" not in df:
            print("{} has no status field".format(f))

        df = df.loc[df["success"] == 0]
        df = df.loc[df["status"] != -1]
        if len(df):
            print(df)

    return

def main():
    prefix = "./simulations/"
    refresh = True

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
        "R10(3900)",
        "R10(4000)",
        "R10(5000)",
        "R10(6000)",
        "R10(7000)",
        "R10(8000)",
    ]
    # directories = [
    #     "R10(1000)",
    #     "R10(2000)",
    #     "R10(3000)",
    #     "R10(4000)",
    #     "R10(6000)",
    #     "R10(7000)",
    #     "R10(8000)",
    #     "LTParameters(4000, Soliton(4000, 3998, 0.999))",
    # ]


    # directories = [
    #     "LTQ{Float64,DT}(4000, Soliton(4000, 3998, 0.9999999701976676))",
    #     "LTQ{Float64,DT}(8000, Soliton(8000, 7998, 0.9999999701976676))",
    # ]

    directories = [
        # "LDPC10(612, 1224, 6810298610506088827)",
        "LDPC10(2400, 4800, 16104329366021199122)",
    ]

    directories = ['R10(1000)', 'RQ(1000)', 'R10(2000)', 'RQ(2000)', 'R10(4000)', 'RQ(4000)']

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
    dfs_success = (load(f, index='overhead', only_success=True) for f in directories)
    dfs_any = (load(f, index='overhead', only_success=False) for f in directories)

    # make plots
    failure_plot(dfs_any, descriptors)
    complexity_plot(dfs_success, descriptors, complexity_multipliers)
    # required_overhead_plot(dfs_any, descriptors)
    # required_complexity_plot(dfs_success, descriptors)
    # cdf_plot(dfs_any, descriptors)
    plt.show()

if __name__ == '__main__':
    main()
