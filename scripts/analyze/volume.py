import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pickle

def read_pickle(dir):
    try:
        with open(dir, 'rb') as file:
            data = pickle.load(file)
            print(f"Contents of {dir}:")
            return data
    except FileNotFoundError:
        print(f"File not found: {dir}")
    except pickle.UnpicklingError:
        print(f"Error unpickling file: {dir}")
    except Exception as e:
        print(f"An error occurred while reading {dir}: {e}")

def find_big_volume(data, dir):
    dic = {}

    # Process data into dic
    for key, value in data.items():
        file_name = key[1].split("/")[-1]

        # Modify file_name if "initial" is in dir and "MSS" is not in the key
        if "initial" in dir and "MSS" not in key[1]:
            file_name_parts = file_name.split("_")
            file_name = "_".join(file_name_parts[0:2]) + "_ERT_" + file_name_parts[-1]
            print(f"Modified file_name: {file_name}")

        mean_volume = value["mean_volume"]
        # Use a set to store unique mean_volume values
        if file_name not in dic:
            dic[file_name] = set()  # Initialize with a set
        dic[file_name].add(mean_volume)

    # Printing the results
    for file_name, volumes in dic.items():
        print(f"File: {file_name}, Mean Volumes: {sorted(volumes)}")

    # Plot dic
    plot_dic(dic, dir)

def plot_dic(dic, dir):
    plt.figure(figsize=(10, 6))

    # Extract data for plotting
    file_names = list(dic.keys())
    mean_volumes = [sorted(list(volumes)) for volumes in dic.values()]

    for i, volumes in enumerate(mean_volumes):
        x_labels = [file_names[i]] * len(volumes)
        plt.scatter(x_labels, volumes, label=file_names[i], alpha=0.7)

    plt.title(f"Mean Volume Distribution - {dir.split('/')[-1]}", fontsize=10)
    plt.xlabel("File Name", fontsize=10)
    plt.ylabel("Mean Volume", fontsize=10)
    plt.xticks(rotation=45, ha="right", fontsize=10)  # Rotate file names for better visibility
    plt.yticks(fontsize=10)  # Adjust y-axis font size
    plt.grid(True)
    plt.tight_layout()

    # Save the plot
    save_path = f"/home/user/SA-EDS/results/{dir.split('/')[-1].replace('.pkl', '.png')}"
    plt.savefig(save_path)
    print(f"Plot saved to {save_path}")
    plt.show()

def get_vol_set(data):
    dic = {}
    temp_lst = ["277", "308", "328", "358"]

    for temp in temp_lst:
        dic[temp] = set()

    for key in data:
        temp = key[0]
        if temp in temp_lst:
            dic[temp].add(data[key]["mean_volume"])

    return dic

if __name__ == '__main__':
    # type_of_l_lst = ["L1", "L2", "L3"]
    type_of_l_lst = ["L1", "L2"]
    iteration_lst = ["initial", "second"]

    # 各データ内の比較をする　特に大きいデータを探す
    for type_of_l in type_of_l_lst:
        for iteration in iteration_lst:
            dir = f"/home/user/SA-EDS/dataset/{type_of_l}_data_int_{iteration}.pkl"
            data = read_pickle(dir)
            find_big_volume(data, dir)

    # データの精度の比較を行う
    for type_of_l in type_of_l_lst:
        # L1_data_int_initial.pklとL1_data_int_second.pklの比較
        plot_data = {}
        for iteration in iteration_lst:
            dir = f"/home/user/SA-EDS/dataset/{type_of_l}_data_int_{iteration}.pkl"
            data = read_pickle(dir)
            vol_data = get_vol_set(data)
            key = dir.split("/")[-1]
            plot_data[key] = vol_data
        print(plot_data)

        data = []

        for file_name, temperatures in plot_data.items():
            for temp, volumes in temperatures.items():
                for volume in volumes:
                    data.append({
                        "Temperature": temp,
                        "Volume": volume,
                        "File": file_name
                    })

        # Create DataFrame from the collected data
        df = pd.DataFrame(data)

        # Create a boxplot
        plt.figure(figsize=(12, 8))
        sns.boxplot(x="Temperature", y="Volume", hue="File", data=df, palette="Set2")

        # Set plot labels and title with fontsize=30
        plt.title("Comparison of Mean Volume Distribution between Initial and Second Data", fontsize=20)
        plt.xlabel("Temperature", fontsize=20)
        plt.ylabel("Mean Volume", fontsize=20)
        plt.xticks(fontsize=20)  # Adjust x-axis font size
        plt.yticks(fontsize=20)  # Adjust y-axis font size
        plt.tight_layout()

        # Show the plot
        plt.show()
        plt.savefig(f"/home/user/SA-EDS/results/{type_of_l}_comp.png")
