#!/usr/bin/env python
"""
@authors: Valerio Vitali, Rebecca Hagen and Francesco Catania (Evolutionary Cell Biology group [IEB]).
SENES: Simulated Evolution of Nuclear Elements Segregation.
Simulates Somatic Assortment of Alleles (e.g.IESs) across asexual generations.
"""


import sys
import os
import argparse
import pysam
import numpy
import scipy.stats
import pandas as pd
from textwrap import dedent
from time import strftime
from subprocess import call
import matplotlib.pyplot as plt
from plotnine import *
from tqdm import tqdm


PYTHON_VERSION = sys.version_info
VERSION = "beta"
PRORAM = "SENES"
AUTHORS = "Valerio Vitali, Rebecca Hagen and Francesco Catania"
CONTACT = "vitaliv@uni-muenster.de"


def core(model, chromosomes, generations, ploidy, allele, input_ratio, out, plot, nullisomics):
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Start modeling of somatic assortment\n")
    sys.stdout.flush()
    #print([model, chromosomes, generations, ploidy, allele, input_ratio, out])
    #day = 10
    #div_24h = 4
    #generations = day * div_24h
    gen_range = range(2, generations + 1, 1)
    if input_ratio and allele is None:
        allele = input_ratio * ploidy
    else:
        input_ratio = allele / ploidy    
    if chromosomes > 1 and input_ratio > 0.5 and nullisomics is False:
        allele = (1 - input_ratio) * ploidy    
    [M, n, N] = [2 * ploidy * chromosomes, 2 * allele, ploidy * chromosomes]
    x = range(0, ploidy * chromosomes + 1, 1) # All possible successes
    y = scipy.stats.hypergeom.pmf(x, M, n, N) # 1st GEN
    df_GEN = []
    
    for j in tqdm(gen_range, desc="progress"):   
        d = [ y[numpy.where(y > 0)[0]][i] * 
              scipy.stats.hypergeom.pmf(x, M, 2 * numpy.where(y > 0)[0][i], ploidy * chromosomes) 
              for i in range(0, len(numpy.where(y > 0)[0]) ) ]        
        y = sum(d) # stack arrays
        if chromosomes > 1 and nullisomics is False:
            x1 = range(0, ploidy + 1)
            y1 = numpy.append(y[0:ploidy], sum(y[ploidy:])) / 2
            #--------------------run magic block
            if input_ratio == 0.5: # sym
                y2 = y1 + y1[::-1]
            elif input_ratio < 0.5:
                y2 = y1*2
            elif input_ratio > 0.5:
                y2 = y1[::-1]*2
            #------------------------------------
        elif chromosomes <= 1 or (chromosomes > 1 and nullisomics is True):
            x1 = x
            y2 = y
        sdev = numpy.sqrt(sum([i**2 * y2[i] for i in x1]) - sum([i for i in x1]*y2)**2)
        if nullisomics is True:
            Hmz = y2[0]
            lim = 3
        else:
            Hmz = y2[0] + y2[ploidy]
            lim = 1
        H = 1 - Hmz # (y2[1:ploidy])
        df_GEN.append(pd.DataFrame({'gen': int(j),
                                    'input_ratio': allele / float(ploidy),
                                    'sdev': sdev,
                                    'H': 1-Hmz,
                                    '1-H': Hmz,
                                    'p(X)': y2}))
    #print(len(x1), len(y2))
    df_GEN = pd.concat(df_GEN)
    df_GEN['X'] = df_GEN.index
    popped_col = df_GEN.pop('X')
    df_GEN.insert(5, 'X', popped_col)
    df_GEN['X'] = df_GEN['X'] / ploidy
    df_GEN['sdev'] = df_GEN['sdev'] / ploidy
    df_GEN['gen'] = df_GEN['gen'].astype('category') # convert gen to factor
    df_GEN_report = df_GEN[ df_GEN['X'] == round(allele / ploidy)].drop_duplicates(subset = 'gen')
    df_GEN_report = df_GEN_report .drop(['X', 'p(X)'], axis=1)
    print(df_GEN_report) # print to stdout
    
    df_GEN.to_csv(path_or_buf = out + "_df_long_" + model[:3] + str(generations), sep = '\t', index = False) # write to file
    df_GEN_report.to_csv(path_or_buf = out + "_df_report_" + model[:3] + str(generations), sep = '\t', index = False)
    
    if plot is True:
        # p(x) dist.
        # group by gen
        p = (ggplot(df_GEN)         # data
         + aes(x='X', y='p(X)', color='gen')    # variables
         + geom_line() # type of plot
         + scale_color_discrete(guide=False)
         + xlim(0, lim)
         + theme_bw(base_size=15)
        )
        
        # sdev change
        df_GEN_report['gen'] = df_GEN_report['gen'].astype('int32') # convert gen to numeric
        df_GEN_report['input_ratio'] = df_GEN_report['input_ratio'].astype('category') # convert input_ratio to categorical
        p1 = (ggplot(df_GEN_report)         # defining what data to use
         + aes(x='gen', y='sdev', color='input_ratio')    # defining what variable to use
         + geom_line() # defining the type of plot to use
         + theme_bw(base_size=15)
        )
        
        # loss of H
        p2 = (ggplot(df_GEN_report)         # defining what data to use
         + aes(x='gen', y='H', color='input_ratio')    # defining what variable to use
         + geom_line() # defining the type of plot to use
         + theme_bw(base_size=15)
        )

        if len(out.split('/')) > 1: # save plot
            ggsave(plot = p, filename = out.split('/')[len(out.split('/'))-1] + "_plot_dist_" + model[:3] + str(generations), path = out.split('/')[0], dpi = 300)
            ggsave(plot = p1, filename = out.split('/')[len(out.split('/'))-1] + "_plot_sdev_" + model[:3] + str(generations), path = out.split('/')[0], dpi = 300)
            ggsave(plot = p2, filename = out.split('/')[len(out.split('/'))-1] + "_plot_H_" + model[:3] + str(generations), path = out.split('/')[0], dpi = 300)
        else:
            ggsave(plot = p, filename = out.split('/')[0] + "_plot_dist_" + model[:3] + str(generations), path = os.getenv("HOME"), dpi = 300)
            ggsave(plot = p1, filename = out.split('/')[0] + "_plot_sdev_" + model[:3] + str(generations), path = os.getenv("HOME"), dpi = 300)
            ggsave(plot = p2, filename = out.split('/')[0] + "_plot_H_" + model[:3] + str(generations), path = os.getenv("HOME"), dpi = 300)

#def cor_rec(model, chromosomes, ploidy, allele, input_ratio, out, num_threads):


def main():
    parser = argparse.ArgumentParser(
        description=dedent('''
        SENES: Simulated Evolution of Nuclear Elements Segregation.
        
        SENES performs a mathematical simulation of Somatic Assortment (SA).
        It currently supports two models of macronuclear architecture.
        
        For a detailed treatment of the mathematical models implemented see:
        Preer, J. R. (1976). Quantitative predictions of random segregation
        models of the ciliate macronucleus. Genetics Research, 27(2), 227-238.
        '''),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-v', '--version', action='version', version='SENES ' + VERSION)
    
    subparsers = parser.add_subparsers(help="You may run SENES in simulation or compare mode.", dest='mode',
                                       description=dedent('''
        SENES runs in two modes.
        For detailed usage of each:
            senes.py module -h
        '''))
    
    parser_sim = subparsers.add_parser('simulator', help="Run with simulation mode") 
    parser_sim.add_argument('-m', '--model', help='Subunit model. One of haploid or chromosomal (default = haploid)', 
                          choices=["haploid", "chromosomal"], default="haploid", required=True)
    
    parser_sim.add_argument('-k', '--ploidy', help='Macronuclear ploidy. Number of genome copies in G1 (after amitosis)',
                          type=int, required=True)
    
    parser_sim.add_argument('-g', '--generations', help='Number of amitotic divisions to simulate',
                          type=int, required=True)    
    
    parser_sim.add_argument('-n', '--allele', help='Initial number of target alleles (A0) in G1 (default = none)',
                            type=int, default = None)
    
    parser_sim.add_argument('-i', '--input_ratio', help='A0 expressed as a fraction of k (Default = None)',
                          type=float, default=None)
    
    parser_sim.add_argument('-c', '--chromosomes', help='Number of somatic chromosomes '
                                                      '(default to 1 for the haploid model)', 
                                                      type=int, default=1)     
    
    parser_sim.add_argument('-o', '--output', help='Output dir (current dir if not specified) and prefix for output file (Default prefix = senes)',
                          default="senes")
    
    parser_sim.add_argument('-p', '--plot', help='Plot distributions and save to file (png). Takes no argument (Default = False)',
                          action="store_true", default=False)
    
    parser_sim.add_argument('--nullisomics', help='All copies of both alleles can be lost (nullisomic loci). Takes no argument (Default = False)',
                          action="store_true", default=False)     
    
    parser_sim.add_argument('-t', '--num_threads', help='Number of threads (Default = 1)', type=int,
                          default=1)
    
    parser_comp = subparsers.add_parser('compare', help="Run with compare mode") 

    args = parser.parse_args()
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.mode == "simulator":
        model = args.model
        chromosomes = args.chromosomes
        ploidy = args.ploidy
        generations = args.generations
        allele = args.allele
        input_ratio = args.input_ratio
        out = args.output
        num_threads = max(args.num_threads, 1)
        plot = args.plot
        nullisomics = args.nullisomics
        
        if model == "chromosomal" and chromosomes == 1:
            print("\nERROR: Please input the number of somatic chromosomes (> 1)\n")
            parser_sim.print_help(sys.stderr)
            sys.exit(1)
            
        if model == "haploid" and chromosomes > 1:
            print("\nWARNING: you selected the haploid model. "
                  "The number of somatic chromosomes will be set to 1\n")
            chromosomes = 1
        
        if model == "haploid" and nullisomics is True:
            print("\nWARNING: you selected the haploid model. "
                  "nullisomics flag will be ignored\n")
            nullisomics = False        
        
        if allele and (allele < 0 or allele > ploidy):
            print("\nERROR: Initial number of target alleles should be between 0 and ploidy\n")
            parser_sim.print_help(sys.stderr)
            sys.exit(1)
            
        if allele is None and input_ratio is None:
            print("\nERROR: Please provide a starting value for either allele or input_ratio\n")
            parser_sim.print_help(sys.stderr)
            sys.exit(1)
            
        if input_ratio and (input_ratio < 0 or input_ratio > 1):
            print("\nERROR: input_ratio should be between 0 and 1\n")
            parser_sim.print_help(sys.stderr)
            sys.exit(1)
            
        if allele and input_ratio:
            print("\nWARNING: Running with allele parameter. Input_ratio will be ignored\n")       

        print("\nRunning SENES with the following shape parameters:\n"
        "--------------------------------------------------")
        print("subunit model:", model)
        if model == "chromosomal":
            print("somatic chromosomes:", chromosomes)
        else: 
            print("somatic chromosomes:", "two haploid whole-genome subunits")
        print("ploidy:", ploidy)
        print("generations:", generations)
        if input_ratio:
            print("allele:", int(round(input_ratio * ploidy)) )
        else:
            print("allele:", allele)
        if allele is None:
            print("input_ratio:", input_ratio )
        else:
            print("input_ratio:", allele / float(ploidy) )
        print("out:", out)
        print("num_threads:", num_threads)
        print("plot:", plot)
        print("nullisomics:", nullisomics)
        print("--------------------------------------------------")        

        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ': ' + ' '.join(sys.argv) + '\n')
        sys.stdout.flush()

        dir_name = os.path.dirname(out)
        if dir_name != '':
            call("mkdir -p " + dir_name, shell=True)

        core(model, chromosomes, generations, ploidy, allele, input_ratio, out, plot, nullisomics)

    elif args.mode == "compare":
        print("\ncomopare mode is under developoment. Watch this space!!\n")
        parser_sim.print_help(sys.stderr)
        sys.exit(1)

    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": DONE!\n")
    sys.stdout.close()


if __name__ == "__main__":
    main()
