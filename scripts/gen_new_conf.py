import common.get_target_file as gtf

# def gen_new_top(top_path, new_top_path):


def main():

    #### 以下を埋める ####  
    # conf_path="/home/user/SA-EDS/sample_structures/L1_r1730290728090-51_277_4/last_conf_r1730290728090-51_277_4.dat.dat"
    # top_path="/home/user/SA-EDS/sample_structures/L1_r1730290728090-51_277_4/generated_r1730290728090-51_277_1.top"
    # req_path="/home/user/SA-EDS/sample_structures/L1_r1730290728090-51_277_4/req_r1730290728090-51_277_4.txt"
    # seq_path="/home/user/SA-EDS/sample_structures/L1_r1730290728090-51_277_4/seq_req_r1730290728090-51_277_4.txt"

    # new_conf_path="/home/user/SA-EDS/sample_structures/L1_r1730290728090-51_277_4/new_last_conf_r1730290728090-51_277_4.dat.dat"
    # new_top_path="/home/user/SA-EDS/sample_structures/L1_r1730290728090-51_277_4/new_generated_r1730290728090-51_277_1.top"

    # complex_lst = [8, 10, 13, 15, 18, 21, 24]
    # seq_length = 20

    # conf_path="/home/user/SA-EDS/sample_structures/L1-GA100000-0.80-ERT-9_277_9/last_conf_L1-GA100000-0.80-ERT-9_277_9.dat.dat"
    # top_path="/home/user/SA-EDS/sample_structures/L1-GA100000-0.80-ERT-9_277_9/generated_L1-GA100000-0.80-ERT-9_277_9.top"
    # seq_path="/home/user/SA-EDS/sample_structures/L1-GA100000-0.80-ERT-9_277_9/seq_req_L1-GA100000-0.80-ERT-9_277_9.txt"

    # new_conf_path="/home/user/SA-EDS/sample_structures/L1-GA100000-0.80-ERT-9_277_9/new_last_conf_L1-GA100000-0.80-ERT-9_277_9.dat.dat"
    # new_top_path="/home/user/SA-EDS/sample_structures/L1-GA100000-0.80-ERT-9_277_9/new_generated_L1-GA100000-0.80-ERT-9_277_9.top"

    # complex_lst = [1, 13, 15, 17, 18, 26, 29, 30]
    # seq_length = 20

    # conf_path="/home/user/SA-EDS/sample_structures/L2-GA100000-0.50-MSS-5_277_7/last_conf_L2-GA100000-0.50-MSS-5_277_7.dat.dat"
    # top_path="/home/user/SA-EDS/sample_structures/L2-GA100000-0.50-MSS-5_277_7/generated_L2-GA100000-0.50-MSS-5_277_7.top"
    # seq_path="/home/user/SA-EDS/sample_structures/L2-GA100000-0.50-MSS-5_277_7/seq_req_L2-GA100000-0.50-MSS-5_277_7.txt"

    # new_conf_path="/home/user/SA-EDS/sample_structures/L2-GA100000-0.50-MSS-5_277_7/new_last_conf_L2-GA100000-0.50-MSS-5_277_7.dat.dat"
    # new_top_path="/home/user/SA-EDS/sample_structures/L2-GA100000-0.50-MSS-5_277_7/new_generated_L2-GA100000-0.50-MSS-5_277_7.top"


    # complex_lst = [2, 3, 7, 11, 14, 18]
    # seq_length = 20

    conf_path="/home/user/SA-EDS/sample_structures/L2_r1730290789727-90_277_7/last_conf_r1730290789727-90_277_7.dat.dat"
    top_path="/home/user/SA-EDS/sample_structures/L2_r1730290789727-90_277_7/generated_r1730290789727-90_277_1.top"
    seq_path="/home/user/SA-EDS/sample_structures/L2_r1730290789727-90_277_7/seq_req_r1730290789727-90_277_7.txt"

    new_conf_path="/home/user/SA-EDS/sample_structures/L2_r1730290789727-90_277_7/new_last_conf_r1730290789727-90_277_7.dat.dat"
    new_top_path="/home/user/SA-EDS/sample_structures/L2_r1730290789727-90_277_7/new_generated_r1730290789727-90_277_1.top"

    complex_lst = [1, 8, 10, 14, 17, 19, 26]
    seq_length = 20

    #####################

    N_lst = []
    for c in complex_lst:
        for i in range(seq_length):
            N_lst.append((c - 1)*seq_length + i + 1)


    seq_f = open(seq_path, "r")

    for i, seq in enumerate(seq_f):
        if i in complex_lst:
            print(i, seq, end="")

    

    conf_f = open(conf_path, "r")
    top_f = open(top_path, "r")
    new_conf_f = open(new_conf_path, "w")
    new_top_f = open(new_top_path, "w")

    # gen new_top

    # The first row contains the total number of nucleotides N and the number of strands Ns:
    N = len(complex_lst)*seq_length
    Ns = len(complex_lst)
    new_top_f.write(f"{N} {Ns}\n")

    num_set = set()

    for i, l in enumerate(top_f):
        if i in N_lst:
            l_lst = l.split(" ")
            str_num = int(l_lst[0])
            num_set.add(str_num)
            new_num = len(num_set)
            l_lst[0] = str(new_num)
            if l_lst[2] != "-1":
                l_lst[2] = str(int(l_lst[2]) - (str_num - 1 - new_num + 1)*seq_length)
            if l_lst[3] != "-1\n":
                l_lst[3] = str(int(l_lst[3]) - (str_num - 1 - new_num + 1)*seq_length) + "\n"
            new_l = " ".join(l_lst)
            new_top_f.write(new_l)

    # gen_new_conf

    # kinetic energies, Etot, U and K, respectively:

    for i, l in enumerate(conf_f):
        if i in [0, 1, 2]:
            new_conf_f.write(l)
        elif i - 2 in N_lst:
            new_conf_f.write(l)


    top_f.close()
    conf_f.close()
    new_top_f.close()
    new_conf_f.close()




if __name__ == '__main__':
    main()