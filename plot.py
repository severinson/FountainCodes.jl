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
    df['overhead'] /= df['K']
    df['num_xor'] /= df['K']
    return df

def failure_plot():
    plt.figure()
    ax = plt.subplot()
    df = load(pattern="*LT*.csv", only_success=False)
    print(df)
    for K in np.unique(df['K']):
        if np.isnan(K):
            continue

        df2 = df[df['K'] == K]
        plt.semilogy(df2['overhead'], df2['success'], 'o', label="LT " + str(K))

    df = load(pattern="*R10*.csv", only_success=False)
    print(df)
    for K in np.unique(df['K']):
        if np.isnan(K):
            continue

        df2 = df[df['K'] == K]
        plt.semilogy(df2['overhead'], df2['success'], 's', label="R10 " + str(K))

    plt.setp(ax.get_xticklabels(), fontsize=25)
    plt.setp(ax.get_yticklabels(), fontsize=25)
    plt.xlabel("Overhead", fontsize=25)
    plt.ylabel("Decoding success probability", fontsize=25)
    plt.xlim([0, 0.5])
    plt.grid()
    plt.legend()
    plt.autoscale(enable=True)
    plt.tight_layout()
    plt.savefig("failure.pdf")

def xor_plot():
    plt.figure()
    ax = plt.subplot()
    ax.set_xlim(xmin=0, xmax=0.5)
    ax.set_ylim(ymin=0, ymax=150)
    df = load(pattern="*LT*.csv")
    print(df)
    # for K in [1000]:
    for K in np.unique(df['K']):
        if np.isnan(K):
            continue

        df2 = df[df['K'] == K]
        plt.plot(df2['overhead'], df2['num_xor'], 'o', label="LT " + str(K))

    df = load(pattern="*R10*.csv")
    print(df)
    # for K in [1000]:
    for K in np.unique(df['K']):
        if np.isnan(K):
            continue

        df2 = df[df['K'] == K]
        plt.plot(df2['overhead'], df2['num_xor'], 's', label="R10 " + str(K))

    plt.setp(ax.get_xticklabels(), fontsize=25)
    plt.setp(ax.get_yticklabels(), fontsize=25)
    plt.xlabel("Overhead", fontsize=25)
    plt.ylabel("XORs/K", fontsize=25)
    plt.grid()
    plt.legend()
    plt.autoscale(enable=True)
    # ax.set_xlim(xmin=0, xmax=0.5)
    ax.set_ylim(ymin=0, ymax=150)
    plt.tight_layout()
    plt.savefig("xor.pdf")

if __name__ == '__main__':
    # failure_plot()
    xor_plot()
    plt.show()
