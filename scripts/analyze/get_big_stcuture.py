import pickle
import sys

def read_pickle(dir):
    try:
        with open(dir, 'rb') as file:
            data = pickle.load(file)
            return data
    except Exception as e:
        print(f"An error occurred while reading the pickle file: {e}")
        return None

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("usage : python3 get_big_structure.py type_of_l target")
        sys.exit(1)

    type_of_l = sys.argv[1]
    target = sys.argv[2]
    base_dir = "/home/user/SA-EDS/"

    dir = f"/home/user/SA-EDS/dataset/{type_of_l}_data_{target}_each.pkl"
    print(dir)
    data = read_pickle(dir)
    if data is not None:
        print("Pickle file loaded successfully.")
        max_val = 0
        max_key = ""
        for key in data:
            try:
                max_nb_strand = data[key]['max_nb_strand']
                print(key, max_nb_strand)
                if max_val < max_nb_strand:
                    max_val = max_nb_strand
                    max_key = key
            except KeyError:
                print(f"KeyError: 'max_nb_strand' not found for {key}")
            except Exception as e:
                print(f"Error in processing {key}: {e}")

        print(max_key, max_val, data[max_key]["max_nb_strand_component"])
    else:
        print("Failed to load the pickle file.")
