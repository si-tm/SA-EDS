import common.get_target_file as gtf

# def gen_new_top(top_path, new_top_path):


def main():

    #### 以下を埋める ####  
    conf_path="/home/user/SA-EDS/sample_structures/L3-GA100000-0.50-ERT-1_277_10/last_conf_L3-GA100000-0.50-ERT-1_277_10.dat.dat"
    top_path="/home/user/SA-EDS/sample_structures/L3-GA100000-0.50-ERT-1_277_10/generated_L3-GA100000-0.50-ERT-1_277_10.top"
    seq_path="/home/user/SA-EDS/sample_structures/L3-GA100000-0.50-ERT-1_277_10/seq_req_L3-GA100000-0.50-ERT-1_277_10.txt"

    new_conf_path="/home/user/SA-EDS/sample_structures/L3-GA100000-0.50-ERT-1_277_10/new_last_conf_L3-GA100000-0.50-ERT-1_277_10.dat.dat"
    new_top_path="/home/user/SA-EDS/sample_structures/L3-GA100000-0.50-ERT-1_277_10/new_generated_L3-GA100000-0.50-ERT-1_277_10.top"

    complex_lst = [8, 17, 20, 22, 28]
    s_lst = []
    base_seq_length_lst = [51, 51, 51, 49, 51, 51]
    

    """
    # structure
    number_of_types = 5
    s0 = a a b @ initial 1.0 M
    s1 = a a c @ initial 1.0 M
    s2 = a a b* @ initial 1.0 M
    s3 = a a d* @ initial 1.0 M
    s4 = a a f* @ initial 1.0 M
    s5 = a b b @ initial 1.0 M

    # length
    length a = 17
    length b = 17
    length c = 17
    length d = 15
    length e = 17
    length f = 17
    """

    #####################


    seq_f = open(seq_path, "r")

    for i, seq in enumerate(seq_f):
        if i + 1 in complex_lst:
            print(i, seq, end="")


    

    conf_f = open(conf_path, "r")
    top_f = open(top_path, "r")
    new_conf_f = open(new_conf_path, "w")
    new_top_f = open(new_top_path, "w")

    # # gen new_top

    # The first row contains the total number of nucleotides N and the number of strands Ns:
    Ns = len(complex_lst)
    for complex in complex_lst:
        if complex%5 == 0:
            s_lst.append(int(complex/5) - 1)
        else:
            s_lst.append(int(complex/5))


    print(complex_lst)
    print(s_lst)

    # 先頭の数
    # base_diff_lst = [357, 814, 961, 1061, 1367]
    base_diff_lst = []

    for i, s in enumerate(s_lst):
        sum = 0
        for j in range(s_lst[i]):
            sum += base_seq_length_lst[j]*5
        if complex_lst[i]%5 != 0:
            base_diff_lst.append((complex_lst[i]%5 - 1)*base_seq_length_lst[s] + sum)
        else:
            base_diff_lst.append(sum + base_seq_length_lst[s_lst[i]]*5 - base_seq_length_lst[s])

    print(base_diff_lst)

    base_diff_lst = [357, 814, 961, 1061, 1367]
    print(base_diff_lst)

    N = 0
    for i, complex in enumerate(complex_lst):
        N += (base_seq_length_lst[s_lst[i]])
    
    new_top_f.write(f"{N} {Ns}\n")

    N_lst = []

    for i, base_diff in enumerate(base_diff_lst):
        # print(base_seq_length_lst[s_lst[i]])
        for j in range(base_seq_length_lst[s_lst[i]]):
            N_lst.append(base_diff + j)
        #     print(base_diff + j, end=" ")
        # print()
       

    s_set = set()

    base_add_lst = []
    base_add_lst.append(0)

    sum = 0

    for s in s_lst:
        sum += base_seq_length_lst[s]
        base_add_lst.append(sum)
    
    for i, l in enumerate(top_f):
        if i - 1 in N_lst:
            l_lst = l.split(" ")
            s_set.add(l_lst[0])
            s_index = len(s_set) - 1
            l_lst[0] = str(s_index + 1)
            if l_lst[2] != "-1":
                l_lst[2] = str(int(l_lst[2]) - int(base_diff_lst[s_index]) + base_add_lst[s_index])
            if l_lst[3] != "-1\n":
                l_lst[3] = str(int(l_lst[3]) - int(base_diff_lst[s_index]) + base_add_lst[s_index]) + "\n"
            new_l = " ".join(l_lst)
            new_top_f.write(new_l)

    # # gen_new_conf

    # # kinetic energies, Etot, U and K, respectively:

    for i, l in enumerate(conf_f):
        if i in [0, 1, 2]:
            new_conf_f.write(l)
        elif i - 3 in N_lst:
            new_conf_f.write(l)


    top_f.close()
    conf_f.close()
    new_top_f.close()
    new_conf_f.close()




if __name__ == '__main__':
    main()