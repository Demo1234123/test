#!/user/bin/env python3
# -*- coding: utf-8 -*-
import re
import os
import pickle
import time
from collections import defaultdict
from typing import DefaultDict, List, Dict, Set
import numpy as np
from aaindex import aaindex1

amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
               'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']


def dynamic_programming_peptides(lengths, property_values, min_value, max_value):
    """
    使用动态规划生成和筛选符合条件的多肽
    :param lengths: 多肽长度范围
    :param property_values: 属性值字典
    :param min_value: 属性值下限
    :param max_value: 属性值上限
    :return: dp和path状态空间
    """
    dp: DefaultDict[int, DefaultDict[int, int]] = defaultdict(lambda: defaultdict(int))
    path: DefaultDict[int, DefaultDict[int, Dict[str, List]]] = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    dp[0][0] = 1  # 初始状态：长度为 0 的多肽，属性值为 0，数量为 1

    # 动态规划计算
    print("开始动态规划计算...")
    start_dp_time = time.time()  # 记录动态规划开始时间
    for length in range(1, max(lengths) + 1):
        for prev_value in dp[length - 1]:
            for aa in amino_acids:
                new_value = prev_value + property_values[aa]
                dp[length][new_value] += dp[length - 1][prev_value]
                # 按氨基酸分类存储路径
                path[length][new_value][aa].append(prev_value)
    end_dp_time = time.time()  # 记录动态规划结束时间
    print(f"动态规划计算完成，用时: {end_dp_time - start_dp_time:.2f} 秒\n")

    return dp, path


def save_single_dp_table(path, index, filename_prefix="output/dp_tables"):
    """
    保存单个动态规划表到本地
    :param path: 单个指标的状态空间
    :param index: 指标索引
    :param filename_prefix: 保存的文件名前缀
    """
    filename = f"{filename_prefix}_{index}.pkl"
    print(f"保存第 {index} 个动态规划表到 {filename}...")
    start_save_time = time.time()

    # 创建输出目录
    os.makedirs(os.path.dirname(filename_prefix), exist_ok=True)

    # 将动态规划表转换为普通字典，以便序列化
    serializable_path = {}
    for length in path:
        serializable_path[length] = {}
        for value in path[length]:
            serializable_path[length][value] = {}
            for aa in path[length][value]:
                serializable_path[length][value][aa] = list(path[length][value][aa])

    # 保存到文件
    with open(filename, 'wb') as f:
        pickle.dump(serializable_path, f)

    end_save_time = time.time()
    print(f"第 {index} 个动态规划表保存完成，用时: {end_save_time - start_save_time:.2f} 秒\n")


def load_all_dp_tables(filename_prefix="output/dp_tables", num_tables=100):
    """
    从本地加载所有动态规划表
    :param filename_prefix: 加载的文件名前缀
    :param num_tables: 要加载的表格数量
    :return: 加载的动态规划表列表和缺失的表格索引列表
    """
    print(f"从 {filename_prefix}_*.pkl 加载动态规划表...")
    start_load_time = time.time()

    criteria_paths = []
    missing_tables = []

    for i in range(num_tables):
        filename = f"{filename_prefix}_{i}.pkl"
        if os.path.exists(filename):
            # 从文件加载
            with open(filename, 'rb') as f:
                serializable_path = pickle.load(f)

            # 将普通字典转换为defaultdict
            path = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
            for length in serializable_path:
                for value in serializable_path[length]:
                    for aa in serializable_path[length][value]:
                        path[int(length)][float(value)][aa] = serializable_path[length][value][aa]

            criteria_paths.append(path)
        else:
            missing_tables.append(i)

    end_load_time = time.time()

    if missing_tables:
        print(f"警告：找不到以下索引的动态规划表: {missing_tables}")
        if not criteria_paths:
            print("没有找到任何动态规划表")
            return None, missing_tables

    print(f"成功加载了 {len(criteria_paths)} 个动态规划表，用时: {end_load_time - start_load_time:.2f} 秒\n")
    return criteria_paths, missing_tables


def find_valid_terminal_states(target_length, criteria_paths, criteria):
    """
    找出目标长度下符合所有指标条件的终端状态和共同氨基酸
    :param target_length: 目标多肽长度
    :param criteria_paths: 所有指标的状态空间
    :param criteria: 指标条件列表
    :return: 有效终端状态字典和共同氨基酸集合
    """
    # 初始化每个指标下的有效氨基酸集合
    valid_aa_sets = []

    # 初始化有效终端状态字典，格式为 {crit_idx: {aa: {curr_value: [prev_value]}}}
    valid_terminal_states = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    # 对于每个指标，找出符合条件的终端状态值和对应的氨基酸
    for i, crit_path in enumerate(criteria_paths):
        min_val = criteria[i]['min_val'] * target_length
        max_val = criteria[i]['max_val'] * target_length

        # 当前指标下的有效氨基酸集合
        valid_aa = set()

        # 找出当前指标下符合条件的终端状态值和对应的氨基酸
        for value in crit_path[target_length]:
            if min_val <= value <= max_val:
                for aa in crit_path[target_length][value]:
                    valid_aa.add(aa)
                    valid_terminal_states[i][aa][value].extend(crit_path[target_length][value][aa])

        valid_aa_sets.append(valid_aa)

    # 计算所有指标下共同的有效氨基酸
    common_aa = set.intersection(*valid_aa_sets) if valid_aa_sets else set()

    # 过滤有效终端状态字典，只保留共同的有效氨基酸
    filtered_valid_terminal_states = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for i in valid_terminal_states:
        for aa in valid_terminal_states[i]:
            if aa in common_aa:
                filtered_valid_terminal_states[i][aa] = valid_terminal_states[i][aa]

    return filtered_valid_terminal_states, common_aa


def backtrack_peptides(target_length, criteria_paths, valid_terminal_states, common_terminal_aa, current_length=1, result=""):
    """
    基于有效终端状态进行回溯
    :param target_length: 目标多肽长度
    :param criteria_paths: 所有指标的状态空间
    :param valid_terminal_states: 有效终端状态字典
    :param common_terminal_aa: 共同的终端氨基酸集合
    :param current_length: 当前回溯的长度
    :param result: 当前构建的多肽序列
    :return: 符合条件的多肽序列列表
    """
    # 如果已经回溯到目标长度，返回结果
    if current_length > target_length:
        return [result]

    peptides = []

    # 如果是第一步回溯（即最后一个氨基酸位置），使用边界条件筛选的共同终端氨基酸
    if current_length == 1:
        for aa in common_terminal_aa:
            # 递归回溯下一个位置
            sub_peptides = backtrack_peptides(target_length, criteria_paths, valid_terminal_states,
                                              common_terminal_aa, current_length + 1, aa)
            peptides.extend(sub_peptides)
    else:
        # 对于中间状态，使用原始的未经筛选的状态空间
        # 获取当前序列的最后一个氨基酸
        last_aa = result[-1]

        # 初始化每个指标下的有效前置氨基酸集合
        valid_prev_aa_sets = []

        # 对于每个指标，找出当前氨基酸对应的前置状态值和前置氨基酸
        for i, crit_path in enumerate(criteria_paths):
            # 当前长度在状态空间中的索引
            curr_length_idx = current_length
            # 前一个长度在状态空间中的索引
            prev_length_idx = current_length - 1

            # 如果是最后一步回溯（即第一个氨基酸位置），需要特殊处理
            if current_length == target_length:
                # 当前氨基酸的所有有效终端状态值
                valid_curr_values = set()
                for value in valid_terminal_states[i][last_aa]:
                    valid_curr_values.update(valid_terminal_states[i][last_aa][value])

                # 当前指标下的有效前置氨基酸集合
                valid_prev_aa = set()

                # 对于当前氨基酸的每个有效终端状态值，找出对应的前置状态值和前置氨基酸
                for prev_value in valid_curr_values:
                    for aa in amino_acids:
                        if prev_length_idx in crit_path and prev_value in crit_path[prev_length_idx]:
                            if aa in crit_path[prev_length_idx][prev_value]:
                                valid_prev_aa.add(aa)

                valid_prev_aa_sets.append(valid_prev_aa)
            else:
                # 当前指标下的有效前置氨基酸集合
                valid_prev_aa = set()

                # 对于当前氨基酸，找出所有可能的前置状态值和前置氨基酸
                for value in crit_path[curr_length_idx]:
                    if last_aa in crit_path[curr_length_idx][value]:
                        for prev_value in crit_path[curr_length_idx][value][last_aa]:
                            for aa in amino_acids:
                                if prev_length_idx in crit_path and prev_value in crit_path[prev_length_idx]:
                                    if aa in crit_path[prev_length_idx][prev_value]:
                                        valid_prev_aa.add(aa)

                valid_prev_aa_sets.append(valid_prev_aa)

        # 计算所有指标下共同的有效前置氨基酸
        common_prev_aa = set.intersection(*valid_prev_aa_sets) if valid_prev_aa_sets else set()

        # 对于每个共同的前置氨基酸，递归回溯
        for prev_aa in common_prev_aa:
            sub_peptides = backtrack_peptides(target_length, criteria_paths, valid_terminal_states,
                                              common_terminal_aa, current_length + 1, prev_aa + result)
            peptides.extend(sub_peptides)

    return peptides


def main(use_cached_dp_tables=True):
    """
    主函数
    :param use_cached_dp_tables: 是否使用缓存的动态规划表
    """
    # 创建输出目录
    os.makedirs("output", exist_ok=True)
    dp_tables_prefix = "output/dp_tables"

    # 加载AAindex数据
    print("加载 AAindex 数据中...")
    start_load_time = time.time()  # 记录加载开始时间
    available_keys = aaindex1.record_codes()[:100]  # 使用前100个AAindex指标
    end_load_time = time.time()  # 记录加载结束时间
    print(f"AAindex 数据加载完成，用时: {end_load_time - start_load_time:.2f} 秒\n")

    # 用户指定的参数
    target_lengths = range(2, 51)
    criteria = []

    # 构建100个测试指标
    for key in available_keys:
        values = aaindex1[key]['values']
        aa_values = [values[aa] for aa in amino_acids]

        # 设置淘汰80%氨基酸的阈值
        lower = np.percentile(aa_values, 10)
        upper = np.percentile(aa_values, 30)

        criteria.append({
            'property': values,
            'min_val': lower,
            'max_val': upper
        })
    # 尝试从本地加载动态规划表
    criteria_paths = []
    missing_tables = list(range(len(criteria)))  # 初始化为所有表格都缺失

    if use_cached_dp_tables:
        criteria_paths, missing_tables = load_all_dp_tables(dp_tables_prefix, len(criteria))

    # 检查是否需要计算缺失的动态规划表
    if missing_tables:
        print(f"需要计算 {len(missing_tables)} 个缺失的动态规划表...")
        start_criteria_time = time.time()

        for i in missing_tables:
            print(f"计算第 {i + 1}/{len(criteria)} 个指标的状态空间...")
            dp, path = dynamic_programming_peptides(target_lengths, criteria[i]['property'], criteria[i]['min_val'], criteria[i]['max_val'])

            # 将新计算的path插入到正确的位置
            if i >= len(criteria_paths):
                # 如果是在列表末尾添加
                criteria_paths.append(path)
            else:
                # 如果是在列表中间插入
                criteria_paths.insert(i, path)

            # 立即保存当前指标的动态规划表
            save_single_dp_table(path, i, dp_tables_prefix)

            # 释放内存
            dp = None

        end_criteria_time = time.time()
        print(f"缺失的状态空间计算完成，总用时: {end_criteria_time - start_criteria_time:.2f} 秒\n")

    # 确保criteria_paths的长度与criteria一致
    if len(criteria_paths) != len(criteria):
        print(f"错误：动态规划表数量 ({len(criteria_paths)}) 与指标数量 ({len(criteria)}) 不匹配！")
        return

    # 对每个目标长度进行处理
    all_filtered_peptides = []

    # 保存筛选结果的文件
    results_file = "output/filtered_peptides.txt"

    for target_length in target_lengths:
        print(f"处理长度为 {target_length} 的多肽...")

        # 计算有效终端状态和共同终端氨基酸
        print(f"计算长度为 {target_length} 的有效终端状态和共同终端氨基酸...")
        valid_terminal_states, common_terminal_aa = find_valid_terminal_states(target_length, criteria_paths, criteria)

        if not common_terminal_aa:
            print(f"长度 {target_length} 没有找到共同终端氨基酸，跳过")
            continue

        print(f"长度 {target_length} 的共同终端氨基酸数量: {len(common_terminal_aa)}")
        print(f"共同终端氨基酸: {', '.join(sorted(common_terminal_aa))}")

        # 回溯生成多肽序列
        print(f"开始回溯生成长度为 {target_length} 的多肽序列...")
        start_backtrack_time = time.time()
        peptides = backtrack_peptides(target_length, criteria_paths, valid_terminal_states, common_terminal_aa)
        end_backtrack_time = time.time()
        print(f"长度 {target_length} 找到的多肽数量: {len(peptides)}")
        print(f"回溯用时: {end_backtrack_time - start_backtrack_time:.2f} 秒")

        all_filtered_peptides.extend(peptides)

        # 将当前长度的结果保存到文件
        with open(f"output/peptides_length_{target_length}.txt", "w") as f:
            for peptide in peptides:
                f.write(f"{peptide}\n")

    # 输出结果
    print(f"筛选完成，符合条件的多肽数量: {len(all_filtered_peptides)}")
    if all_filtered_peptides:
        print("示例符合条件的多肽:")
        print(all_filtered_peptides[:min(10, len(all_filtered_peptides))])  # 打印前10个结果

        # 保存所有筛选结果到文件
        with open(results_file, "w") as f:
            for peptide in all_filtered_peptides:
                f.write(f"{peptide}\n")
        print(f"所有筛选结果已保存到 {results_file}")
    else:
        print("没有找到符合所有条件的多肽。")


if __name__ == "__main__":
    main(use_cached_dp_tables=True)  # 默认使用缓存的动态规划表
