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
    ave_stability = 0
    num = 0

    if data is not None:
        print("Pickle file loaded successfully.")
        max_val = 0
        max_key = ""
        for key in data:
            if "sigmoid" in data[key] and key[1] == "277":
                try:
                    print(key[1],data[key]["sigmoid"]["b"], data[key]["sigmoid"]["c"][2])
                    b = float(data[key]["sigmoid"]["b"])
                    c = float(data[key]["sigmoid"]["c"][2])
                    if c != 0:
                        ave_stability += b/c
                        num += 1
                except:
                    print(key[1],data[key]["sigmoid"]["b"], data[key]["sigmoid"]["c"])
                    b = float(data[key]["sigmoid"]["b"])
                    c = float(data[key]["sigmoid"]["c"])
                    if c != 0:
                        ave_stability += b/c
                        num += 1
    else:
        print("Failed to load the pickle file.")

    ave_stability /= num
    print(ave_stability)
