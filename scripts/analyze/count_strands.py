import sys
sys.path.append("common/")
sys.path.append("measuring_volume/")
import measuring_volume.run_output_bonds as rob
import common.get_target_file as gtf
# import convexhull_volume as cv
import measuring_volume.get_top_data as gtd

# strand数の前後を書く．
def count_strands(target_dir):
    strands2particle, particle2strand = gtd.make_initial_strands_data(target_dir)
    # print("before : ", len(strands2particle))
    before = len(strands2particle)
    strands2particle, particle2strand = gtd.get_connection_strands(gtf.get_bonds(target_dir), strands2particle, particle2strand)
    # print("after : ", len(strands2particle))
    after = len(strands2particle)

    return (before, after)

# strand数の前後を書く．
def analyze_strands(target_dir):
    strands2particle, particle2strand = gtd.make_initial_strands_data(target_dir)
    strands2particle, particle2strand = gtd.get_connection_strands(gtf.get_bonds(target_dir), strands2particle, particle2strand)

    # which paricle are bonded
    # for strand in strands2particle:
        # print(strand)
        # print(strand, strands2particle[strand])
    # print(particle2strand)

# 最大結合strand数を抽出する
def max_nb_strands(target_dir):
    bf_strands2particle, bf_particle2strand = gtd.make_initial_strands_data(target_dir)
    strands2particle = bf_strands2particle
    particle2strand = bf_particle2strand
    str_key = print("number ", bf_strands2particle.keys())

    af_strands2particle, af_particle2strand = gtd.get_connection_strands(gtf.get_bonds(target_dir), strands2particle, particle2strand)

    print("number ", bf_strands2particle.keys())
    # which paricle are bonded
    for strand in bf_strands2particle:
        print(strand, len(bf_strands2particle[strand]))
        if strand in af_strands2particle:
            print(strand, len(af_strands2particle[strand]))
        else:
            print("none")
    # print(strands2particle)


def test():
    target_dir = "../input/results/initial/L3_initial_6/L3-GA100000-0.50-ERT-6_277_6"
    before, after = count_strands(target_dir)
    # print(before, after)
    # analyze_strands(target_dir)
    max_nb_strands(target_dir)
    print(after/before)

# def main():

if __name__ == '__main__':
    # main()
    test()