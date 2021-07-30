import os


def get_list(file):
    """list of all data files"""
    samples_lst = []
    with open(file) as fi:
        for line in fi:
            samples_lst.append(line.strip())
    return samples_lst


def get_sample_ids(samples_lst):
    """list of just sample prefix"""
    sample_ids = []
    for sample in samples_lst:
        sample_ids.append((sample.split("_")[0]))
    sorted_sample_ids = sorted(set(sample_ids), key=lambda x: float("." + x[3:]))
    print(sorted_sample_ids)
    print(len(sorted_sample_ids) == len(set(sorted_sample_ids)))
    return sorted_sample_ids


def get_sample_sequence_ids(samples_lst):
    """list of just sample prefix"""
    sample_sequence_ids = []
    for sample in samples_lst:
        sample_sequence_ids.append(("_".join(sample.split("_")[0:2])))
    sorted_sample_sequence_ids = sorted(set(sample_sequence_ids), key=lambda x: float("." + x[3:8]))
    return sorted_sample_sequence_ids


def make_samples_dict(sample_ids):
    """helper function for get_sample_pairs() - dict of sample id with r1/r2 value pairs as [str, str]"""
    samples_dict = {}
    for sample_id in sample_ids:
        samples_dict[sample_id] = ["", ""]
    return samples_dict


def get_sample_pairs(samples_lst, sample_ids):
    """dict of sample id with [r1, r2] as value pair"""
    samples_dict = make_samples_dict(sample_ids)
    for sample_id in sample_ids:
        cnt = 0
        for sample in samples_lst:
            if sample_id in sample:
                if "R1" in sample:
                    samples_dict[sample_id][0] = sample
                    cnt += 1
                elif "R2" in sample:
                    samples_dict[sample_id][1] = sample
                    cnt += 1
            if cnt == 2:
                break
            elif cnt > 2:
                print("There's an error in the code")
    if len(sample_ids) == len(samples_dict):
        print("Dope, good to go")
    return samples_dict


def merge_dicts(samples_dict_1, samples_dict_2):
    """merge dicts of pairs so {sample: [1.r1, 1.r2, 2.r1, 2.r2]"""
    samples_dict_new = {}
    for key in samples_dict_1:
        samples_dict_new[key] = samples_dict_1[key] + samples_dict_2[key]
    return samples_dict_new


"""
def write_list(samples_dict_new, path, wildcard_name):
    with open("../../config/config.yml", "a+") as fi:
        fi.write("samples:\n")
        fi.write("")
        for key in samples_dict_new:
            fi.write("    {}:\n".format(key))
            fi.write("        {}:\n".format("first"))
            fi.write("            - {}\n".format(samples_dict_new[key][0]))

            fi.write("            - {}\n".format(samples_dict_new[key][1]))

            fi.write("        {}:\n".format("second"))
            fi.write("            - {}\n".format(samples_dict_new[key][2]))

            fi.write("            - {}\n".format(samples_dict_new[key][3]))
"""


def write_samples_tsv_merged(sample_ids):
    """makes a samples tsv file with staggered reads"""
    # this is the bullshit line
    samples_ids_changed = []
    for i in sample_ids:
        samples_ids_changed.append("_".join(i.split("_")[0:2]))
    sample_ids = set(samples_ids_changed)
    with open("../../config/samples_merged.tsv", "w+") as fi:
        fi.write("sample\n")
        for i in sample_ids:
            fi.write("{}\n".format(i))


def write_samples_tsv(sample_ids, suffix):
    """make sample tsv file of samples, specific to the round"""
    # this is the bullshit line
    samples_ids_changed = []
    for i in sample_ids:
        samples_ids_changed.append("_".join(i.split("_")[0:2]))
    sample_ids = set(samples_ids_changed)
    tmp = "../../config/{}samples.tsv".format(suffix)
    with open(tmp, "w+") as fi:
        fi.write("sample\n")
        for i in sample_ids:
            fi.write("{}\n".format(i))


def write_units(samples_dict_new):
    with open("../../config/units.tsv", "w+") as fi:
        fi.write("sample\tunit\tfq1\tfq2\n")
        for key in samples_dict_new:
            fi.write("{}\t1\t{}\t{}\n".format(key, samples_dict_new[key][0], samples_dict_new[key][1]))
            fi.write("{}\t2\t{}\t{}\n".format(key, samples_dict_new[key][2], samples_dict_new[key][3]))


def write_samples_v2_merged(samples_dict_new):
    with open("../../config/samples2.tsv", "w+") as fi:
        fi.write("sample\tf_fq1\tf_fq2\ts_fq1\ts_fq2\n")
        for key in samples_dict_new:
            fi.write("{}\t{}\t{}\t{}\t{}\n".format(key,
                                                   samples_dict_new[key][0], samples_dict_new[key][1],
                                                   samples_dict_new[key][2], samples_dict_new[key][3]))


def main():

    # Write tsv of sample ids with sequencing ids from round 1
    fname_1 = "../../resources/first_round_read_IDs.txt"
    samples_lst_1 = get_list(file=fname_1)
    sample_ids_1 = get_sample_ids(samples_lst_1)
    samples_dict_1 = get_sample_pairs(samples_lst_1, sample_ids_1)

    # Write tsv of only sample names
    write_samples_tsv(sample_ids_1, "")

    sample_sequence_ids_1 = get_sample_sequence_ids(samples_lst_1)
    write_samples_tsv(sample_sequence_ids_1, "first_")

    # Write tsv of sample ids with sequencing IDs (i.e., NTM#####_S##)  from round 2
    fname_2 = "../../resources/second_round_read_IDs.txt"
    samples_lst_2 = get_list(file=fname_2)
    sample_ids_2 = get_sample_ids(samples_lst_2)
    samples_dict_2 = get_sample_pairs(samples_lst_2, sample_ids_2)
    sample_sequence_ids_2 = get_sample_sequence_ids(samples_lst_2)
    write_samples_tsv(sample_sequence_ids_2, "second_")

    # samples_lst_new = samples_lst_1 + samples_lst_2
    samples_dict_new = merge_dicts(samples_dict_1, samples_dict_2)
    # write_samples_tsv_merged(samples_lst_new)
    write_units(samples_dict_new)
    # write_list(samples_dict_new, path_2, wildcard_name_2)

    # write samples v2 merged
    write_samples_v2_merged(samples_dict_new)


if __name__ == '__main__':
    main()
