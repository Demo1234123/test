#!/user/bin/env python3
# -*- coding: utf-8 -*-
import os
import glob
import tarfile


def compress_faa_files(directory, output_file):
    faa_ls = glob.glob(os.path.join(directory, '**/*.faa'), recursive=True)
    with tarfile.open(output_file, "w:gz") as tar:
        for faa_pa in faa_ls:
            if os.path.getsize(faa_pa) <= 2 * 1024 * 1024 * 1024:
                tar.add(faa_pa, arcname=os.path.relpath(faa_pa, directory))
    print(f"所有小于2G的faa文件已打包到 {output_file}")


if __name__ == '__main__':
    base_dir = 'down'
    compress_faa_files(base_dir, 'less_2G_protein.tar.gz')
