#!/user/bin/env python3
# -*- coding: utf-8 -*-
import os
import tarfile
import gzip
import shutil
from Bio import SeqIO
from tqdm import tqdm


def tar_gz(ls):
    # 遍历每一个 .tar 文件
    for pro in tqdm(ls):
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
                        get_seq(extract_path)


def get_seq(seq_pa):
    sequences = list(SeqIO.parse(seq_pa, "fasta"))
    # 筛选长度小于等于50的序列
    filtered_sequences = []
    for seq_record in sequences:
        seq_length = len(seq_record.seq)
        if seq_length <= 50:
            filtered_sequences.append(seq_record)
        elif seq_record.seq.startswith("*") or seq_record.seq.endswith("*"):
            clean_seq_length = len(seq_record.seq.replace('*', ''))
            if clean_seq_length <= 50:
                filtered_sequences.append(seq_record)
    with open(seq_pa, "w") as output_handle:
        SeqIO.write(filtered_sequences, output_handle, "fasta")


if __name__ == '__main__':
    # 目标目录
    base_dir = 'down'
    # 列出目录中的所有文件
    protein_ls = os.listdir(base_dir)
    tar_gz(protein_ls)
