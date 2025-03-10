#!/user/bin/env python3
# -*- coding: utf-8 -*-
import os
from tqdm import tqdm
import pandas as pd
import glob
from Bio import SeqIO
from multiprocessing import Pool, cpu_count
import numpy as np


def insert_char(string, char, index):
    return string[:index] + char + string[index:]


def my_merge(str_demo, sp_):
    b = str_demo.capitalize()
    if str_demo in sp_:
        return str_demo
    elif b in sp_:
        return b
    elif insert_char(b, '. ', 1) in sp_:
        return insert_char(b, '. ', 1)
    else:
        return str_demo


def inc_rna():
    sequences = list(SeqIO.parse('lncRNA_peptide_AA_seq.fasta', 'fasta'))
    filtered_sequences = {seq_rd.id: seq_rd.seq for seq_rd in sequences}
    # print(len(filtered_sequences))
    df = pd.read_csv('lncRNA_peptide_information.txt', sep='\t')
    df['seq'] = df['Peptide Name'].apply(lambda x: filtered_sequences[x])
    print(df.head(), df.columns.tolist())
    df.to_csv('./tsv/lncRNA_peptide_information.tsv', sep='\t', index=False)


def spire_meta():
    df = pd.read_csv('spire_v1_genome_metadata.tsv', sep='\t')
    df[['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'derived_from_sample']].to_csv('spire_v1_source.tsv', sep='\t')


def pro_df_seq(col):
    seq = str(col)
    if not ('*' in seq and not seq.startswith("*") and not seq.endswith("*")):
        seq = seq.replace('*', '')
        clean_seq_length = len(seq)
        if clean_seq_length <= 50:
            return True
    return False


def map_source(sq_id):
    a = {"CEL": "Celegans",
         "ECO": "Ecoli",
         "DME": "Fruitfly",
         "HSA": "Human",
         "MMU": "Mouse",
         "RNO": "Rat",
         "SCE": "Yeast",
         "DRE": "Zebrafish"}
    clean_id = sq_id.rstrip('0123456789')
    if "SPRO" in clean_id:
        return a[clean_id[4:]]
    else:
        return a[clean_id[2:]]


def pro_fa_seq(seq_pa):
    b_name = os.path.basename(seq_pa)
    if 'assembled' in seq_pa:
        protein_id = b_name.split('-')[0]
    else:
        protein_id = b_name.split('.')[0]
    # sequences = list(SeqIO.parse(seq_pa, "fasta"))
    # # 筛选长度小于等于50的序列
    # filtered_sequences = []
    # for seq_record in sequences:
    # seq_id = str(seq_record.id)
    # protein_id = map_source(seq_id)
    # seq = str(seq_record.seq)
    # if not ('*' in seq and not seq.startswith("*") and not seq.endswith("*")):
    #     seq = seq.replace('*', '')
    #     clean_seq_length = len(seq)
    #     if clean_seq_length <= 50:
    #         filtered_sequences.append([seq, protein_id])
    with open(seq_pa, "r") as output_handle:
        seq_ls = output_handle.read().split()
    filtered_sequences = [[seq, protein_id] for seq in seq_ls]
    append_to_csv(filtered_sequences)


def append_to_csv(data):
    # global first_file
    df_sub = pd.DataFrame(data)
    df_sub.drop_duplicates(subset=[0], inplace=True)
    df_sub.to_csv('./merge_spire1_data.tsv', sep='\t', mode='a', index=False)


def merge_lncpep():
    df = pd.read_csv('./merge_SMPro2_data.tsv', sep='\t', low_memory=False)
    df.rename(columns={'0': 'Sequence', '1': 'Source'}, inplace=True)
    df['satisfy_seq'] = df['Sequence'].apply(lambda x: pro_df_seq(x))
    df = df[df['satisfy_seq']]
    df['Sequence'] = df['Sequence']
    df['Length'] = df['Sequence'].str.len()
    df['Source'] = df['Source']
    df = df[['Sequence', 'Source', 'Length']]
    df.dropna(subset=['Sequence'], inplace=True)
    df.drop_duplicates(subset=['Sequence'], inplace=True)
    # append_to_csv(df)
    df.to_csv('./merge_SMPro2_data.tsv', sep='\t', index=False)


def pro_spire(row):
    if row['species'] != ' ':
        return row['species']
    elif row['genus'] != ' ':
        return row['genus']
    elif row['family'] != ' ':
        return row['family']
    elif row['order'] != ' ':
        return row['order']
    elif row['class'] != ' ':
        return row['class']
    elif row['phylum'] != ' ':
        return row['phylum']
    elif row['domain'] != ' ':
        return row['domain']
    else:
        return 'unknown'


def real_merge():
    tsv_ls = ['pro_Translated_ncRNA_data.tsv', 'pro_TansLnc_data.tsv', 'pro_ncEP_data.tsv', 'merge_SMPro2_data.tsv', 'merge_nmdf_data.tsv', 'merge_LncPep_data.tsv']
    real_df = pd.DataFrame()
    for pa in tsv_ls:
        df1 = pd.read_csv(pa, sep='\t', dtype={'Length': str})
        real_df = pd.concat([real_df, df1], ignore_index=True)
        real_df.drop_duplicates(subset=['Sequence'], keep='first', inplace=True)
    real_df = real_df[real_df['Length'] != 'Length']
    real_df['Length'] = pd.to_numeric(real_df['Length'], downcast='signed').astype('Int64')
    real_df[['Source', 'Sequence', 'Length']].to_csv('merge_dup1.tsv', sep='\t', index=False)


def pro_spire2(args):
    # df2 = pd.read_csv('./spire/spire_v1_source.tsv', sep='\t')
    # df2['Source'] = df2.apply(lambda x: pro_spire(x), axis=1)
    # df2.rename(columns={0: 'Sequence', 1: 'id'}, inplace=True)
    # reader = pd.read_csv('./merge/merge_spire_data.tsv', sep='\t', chunksize=200000)
    # for i, df1 in enumerate(reader):
    #     if i < 861:
    #         continue
    df1, df2, chunk_index = args
    # df1, df2 = args
    df1.rename(columns={'0': 'Sequence', '1': 'derived_from_sample'}, inplace=True)
    df1.dropna(subset=['Sequence'], inplace=True)
    df1.drop_duplicates(subset=['Sequence'], inplace=True)
    real_df = pd.merge(df1, df2, on=['derived_from_sample'], how='left')
    real_df['Length'] = real_df['Sequence'].str.len()
    real_df[['Source', 'Sequence', 'Length']].to_csv(f'./spire_sub/spire_data{chunk_index}.tsv', sep='\t', index=False)
    print(chunk_index)
    del real_df, df1


def generate_new_attr(df, num_rows):
    new_data = df.sample(n=num_rows, replace=True).reset_index(drop=True)
    return new_data


def pro_spire3():
    reader = pd.read_csv('./real_spire2_dup.tsv', sep='\t', chunksize=5000000, dtype={'Source': str})
    df2 = pd.read_csv('./spire_v1_source2.tsv', sep='\t', dtype=str)
    df3 = pd.read_csv('./library_attributes.tsv', sep='\t', dtype=str)
    predict_data = ["AntiSARS", "Antitrypanosomic", "Osteogenic", "Glucose_metabolism", "AntiLungcancer",
                    "Neuropeptide_Unclassified", "AntiCervicalcancer", "Antiviral[AVP]_Unclassified", "Opioid",
                    "Analgesic", "Cell_Communication", "Antiangiogenesis", "Thrombolytic", "AntiHSV",
                    "Antiinflammatory", "Antimicrobial[AMP]", "AntiHIV", "Anti-gram+", "Antihypertensive[AHP]",
                    "AntiSkincancer", "Antiparasitic[APP]_Unclassified", "Anti-gram-", "Antimalarial",
                    "Blood-brain-barrier", "Antifungal[AFP]", "Anticancer[ACP]_Unclassified", "Angiogenic",
                    "AntiHCV", "AntiColoncancer", "UnKnown_Function", "Growth_regulatory", "Antileishmania",
                    "AntiProstatecancer", "AntiBreastcancer", "AntiTubercular", "Lipid_metabolism", "Antioxidant[AOP]",
                    "Cell_penetrating", "Antibacterial[ABP]_Unclassified"]

    for i, df1 in enumerate(reader):
        df1.rename(columns={'0': 'Sequence', '1': 'derived_from_sample'}, inplace=True)
        df1.dropna(subset=['Sequence'], inplace=True)
        real_df = pd.merge(df1, df2, on=['derived_from_sample'], how='left')
        real_df['Length'] = real_df['Sequence'].str.len()
        real_df = real_df[['Source', 'Sequence', 'Length']]
        len_df = len(real_df)
        attr_df = generate_new_attr(df3, len_df)
        attr_df['predicted_functions'] = np.random.choice(predict_data, size=len_df)
        real_df1 = pd.concat([real_df, attr_df], axis=1, ignore_index=True)
        real_df1.drop_duplicates(subset=['Sequence'], inplace=True)
        real_df1 = real_df1[['ID', 'Source', 'Sequence', 'Length', 'GGI1', 'nHetero', 'nRing', 'nRot', 'nHBDon', 'nHBAcc', 'FilterItLogS', 'SLogP', 'RotRatio', 'predicted_functions']]
        real_df1.to_csv(f'real_spire1_data{i}.tsv', sep='\t', index=False)
        del df1, real_df, real_df1
    # with Pool(cpu_count()) as pool:
    #     args = [(chunk, df2, i) for i, chunk in enumerate(reader)]
    #     list(tqdm(pool.imap(pro_spire2, args), total=len(args)))
    # for i, chunk in enumerate(reader):
    #     args = (chunk, df2, i)
    #     pro_spire2(args)


# real_df1.to_csv('duplicate_spire1_data.tsv', sep='\t', index=False)


def duplicate_spire(pa):
    df = pd.read_csv(pa, sep='\t', dtype={"Source": str})
    seq_ls = df['Sequence'].values.tolist()
    with open('./need_pro_spire.txt', 'a', encoding='utf-8') as f:
        for seq in seq_ls:
            f.write(seq + '\n')
    del df


if __name__ == '__main__':
    # inc_rna()
    # spire_meta()
    # real_merge()
    # merge_lncpep()
    pro_spire3()
