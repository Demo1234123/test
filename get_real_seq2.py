#!/user/bin/env python3
# -*- coding: utf-8 -*-
import glob
import os
import tarfile
import gzip
import shutil
from Bio import SeqIO
import pandas as pd
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
import subprocess


def process_tar(args):
    pro, base_dir = args
    # 遍历每一个 .tar 文件
    # for pro in ls:
    if pro.endswith('.tar'):
        tar_path = os.path.join(base_dir, pro)

        # 创建一个目录用于解压 .tar 文件
        extract_dir = os.path.join(base_dir, pro[:-4])
        os.makedirs(extract_dir, exist_ok=True)

        # 解压 .tar 文件
        with tarfile.open(tar_path, "r:") as tar:
            tar.extractall(path=extract_dir)

        # 遍历解压后的目录，查找 .gz 文件并解压
        for root, dirs, files in os.walk(extract_dir):
            for file in files:
                if file.endswith('.gz'):
                    gz_path = os.path.join(root, file)
                    extract_path = os.path.join(root, file[:-3])  # 去掉 .gz 后缀
                    # 解压 .gz 文件
                    with gzip.open(gz_path, 'rb') as f_in:
                        with open(extract_path, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    # 可选择删除原始的 .gz 文件
                    os.remove(gz_path)
                    out_path = extract_path + 'a'
                    # get_seq(extract_path)
                    subprocess.run(['prodigal', '-i', extract_path, '-a', out_path])
                    os.remove(extract_path)


def get_seq(seq_pa):
    sequences = list(SeqIO.parse(seq_pa, "fasta"))
    # 筛选长度小于等于50的序列
    filtered_sequences = []
    for seq_record in sequences:
        seq = str(seq_record.seq)
        # seq_length = len(seq)
        # if seq_length <= 50 and '*' not in seq:
        #     filtered_sequences.append(seq)
        if not ('*' in seq and not seq.startswith("*") and not seq.endswith("*")):
            seq = seq.replace('*', '')
            clean_seq_length = len(seq)
            if clean_seq_length <= 50:
                filtered_sequences.append(seq)
    with open(seq_pa, "w") as output_handle:
        for seq in filtered_sequences:
            output_handle.write(seq + "\n")


def filter_sequences(sequence):
    seq_length = len(sequence)
    if seq_length <= 50:
        return True
    elif sequence.startswith("*") or sequence.endswith("*"):
        clean_seq_length = len(sequence.replace('*', ''))
        if clean_seq_length <= 50:
            return True
    return False


# def pro_txt(pa):
#     df = pd.read_csv(pa, sep='\t')
#     if 'orf_sequence' in df.columns.tolist():
#         df['orf_sequence'] = df['orf_sequence'].astype(str)
#         df['is_filtered'] = df['orf_sequence'].apply(filter_sequences)
#     else:
#         df['Protein or peptide sequence'] = df['Protein or peptide sequence'].astype(str)
#         df['is_filtered'] = df['Protein or peptide sequence'].apply(filter_sequences)
#     filtered_df = df[df['is_filtered']]
#     b_name = os.path.basename(pa)
#     filtered_df.to_csv(f'./RNA/pro/{b_name}', sep='\t', index=False)


if __name__ == '__main__':
    # 目标目录
    base_dir = 'test'
    # 列出目录中的所有文件
    protein_ls = os.listdir(base_dir)
    # os.makedirs(os.path.join(base_dir, 'pro'), exist_ok=True)
    # protein_ls = glob.glob(os.path.join(base_dir, '**/*.txt'), recursive=True)
    # pro_txt('pep.info.ZMA.txt')
    # for pa in protein_ls:
    #     get_seq(pa)
    #     pro_txt(pa)
    #     print(pa)
    with Pool(cpu_count()) as pool:
        args = [(pro, base_dir) for pro in protein_ls]
        list(tqdm(pool.imap(process_tar, args), total=len(protein_ls)))
