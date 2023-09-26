import argparse

def main(args):
    import pandas as pd
    import scipy.stats
    import numpy as np

    df1_list = [pd.read_csv(i,sep = "\t",header=0) for i in args.files1]
    dict1 = dict(zip(args.files1,df1_list))

    df2_list = [pd.read_csv(i,sep = "\t",header=0) for i in args.files2]
    dict2 = dict(zip(args.files2,df2_list))

    gene_list = []
    JS_list = []
    for gene in args.genes:
        tmp_JS_list = []
        for file1 in args.files1:
            if gene not in dict1[file1].columns:
                print(f'{gene} not in {file1}')
                continue
            prob1 = np.array(dict1[file1][gene])
            if np.sum(prob1)==0:
                continue
            prob1 = prob1 / np.sum(prob1)
            for file2 in args.files2:
                if gene not in dict2[file2].columns:
                    print(f'{gene} not in {file2}')
                    continue
                prob2 = np.array(dict2[file2][gene])
                if np.sum(prob2) == 0:
                    continue
                prob2 = prob2/np.sum(prob2)
                M = (prob1 + prob2) / 2
                JS_result = 0.5 * scipy.stats.entropy(prob1, M, base=2) + 0.5 * scipy.stats.entropy(prob2, M, base=2)
                tmp_JS_list = tmp_JS_list + [JS_result]
        if len(tmp_JS_list)==0:
            print(f'all files has not {gene}')
            continue
        gene_list = gene_list + [gene]
        JS_list = JS_list + [np.mean(tmp_JS_list)]
    save = pd.DataFrame({'gene':gene_list,args.type:JS_list})
    save.to_csv(f'{args.output}.csv',index=False,header=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f1", "--files1", help="file1 file1 ...", type=str, default=[],nargs='+')
    parser.add_argument("-f2", "--files2", help="file2 file2 ...", type=str, default=[],nargs='+')
    parser.add_argument("-g", "--genes", help="-g gene1 gene2 ...", type=str, default=[], nargs='+')
    parser.add_argument("-t", "--type", help="type  will be the colname ...", type=str, default='JS')
    parser.add_argument("-o", "--output", help="output img name, default output", type=str, default='output')
    args = parser.parse_args()
    main(args)



