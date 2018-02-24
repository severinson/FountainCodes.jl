import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from os import path

plt.rc('pgf',  texsystem='pdflatex')
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{lmodern}']

def success_filter(df):
    return df.loc[df['success'] == 1.0, :]

def mvc(K):
    return K*(K-1)

def load(pattern='*LT*.csv', only_success=True, directory='results'):
    '''load df's from directory, average the columns and return a new df'''
    pattern = path.join(directory, pattern)
    results = list()
    for df in [pd.read_csv(f) for f in glob.glob(pattern)]:
        if only_success:
            df = success_filter(df)
        result = dict()
        for column in df:
            result[column] = df[column].mean()
        results.append(result)

    df = pd.DataFrame(results)
    if len(df):
        df['overhead'] /= df['K']
        df['num_xor'] /= df['K']
        # df['num_xor'] /= mvc(df['K'])
    return df

def failure_plot():
    plt.figure()
    ax = plt.subplot()
    df = load(pattern="*LT*.csv", only_success=False)
    if len(df):
        print(df.head())
        for K in np.unique(df['K']):
            if np.isnan(K):
                continue
            K = int(K)

            df2 = df[df['K'] == K]
            plt.semilogy(df2['overhead'], df2['success'], 'o', label="LT, " + str(K))

    df = load(pattern="*R10*.csv", only_success=False)
    if len(df):
        print(df.head())
        for K in np.unique(df['K']):
            if np.isnan(K):
                continue
            K = int(K)

            df2 = df[df['K'] == K]
            plt.semilogy(df2['overhead'], df2['success'], 's', label="R10, " + str(K))

    plt.setp(ax.get_xticklabels(), fontsize=25)
    plt.setp(ax.get_yticklabels(), fontsize=25)
    plt.xlabel("Overhead", fontsize=25)
    plt.ylabel("Decoding success probability", fontsize=25)
    # plt.xlim([0, 0.5])
    plt.grid()
    plt.legend(
        fontsize=18,
        ncol=1,
        labelspacing=0,
        columnspacing=0.05,
        borderaxespad=0.1,
    )
    plt.autoscale(enable=True)
    plt.tight_layout()
    plt.savefig("failure.pdf")
    plt.savefig("failure.png")

def xor_plot():
    plt.figure()
    ax = plt.subplot()
    ax.set_xlim(xmin=0, xmax=0.5)
    ax.set_ylim(ymin=0, ymax=150)
    df = load(pattern="*LT*.csv")
    if len(df):
        print(df.head())
        for K in np.unique(df['K']):
            if np.isnan(K):
                continue
            K = int(K)

            df2 = df[df['K'] == K]
            plt.plot(df2['overhead'], df2['num_xor'], 'o', label="LT, " + str(K))

    df = load(pattern="*R10*.csv")
    if len(df):
        print(df.head())
        for K in np.unique(df['K']):
            if np.isnan(K):
                continue
            K = int(K)

            df2 = df[df['K'] == K]
            plt.plot(df2['overhead'], df2['num_xor'], 's', label="R10, {}".format(K))

    plt.setp(ax.get_xticklabels(), fontsize=22)
    plt.setp(ax.get_yticklabels(), fontsize=22)
    plt.xlabel(r"$\rm{Overhead}$", fontsize=22)
    # plt.ylabel(r"$\rm{XORs}/(m(m-1))$", fontsize=22)
    plt.ylabel(r"$\rm{XORs}/m$", fontsize=22)
    plt.grid()
    plt.legend(
        fontsize=18,
        ncol=2,
        labelspacing=0,
        columnspacing=0.05,
        borderaxespad=0.1,
        loc="upper left"
    )
    plt.autoscale(enable=True)
    ax.set_xlim(xmin=0.1, xmax=0.5)
    # ax.set_ylim(ymin=20, ymax=160)
    # ax.set_xlim(xmin=0.2)
    ax.set_ylim(ymin=0, ymax=40)
    plt.tight_layout()
    plt.savefig("xor_opt.pdf")
    plt.savefig("xor_opt.png")

if __name__ == '__main__':
    failure_plot()
    xor_plot()
    plt.show()
