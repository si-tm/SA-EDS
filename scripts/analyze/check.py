import pickle
import matplotlib.pyplot as plt

def plot_and_save(dic, output_path):
    """
    データをプロットして保存する
    """
    plt.figure(figsize=(10, 6))
    for key, values in dic.items():
        plt.plot(values, label=f"L1 {key}", marker='o')
    
    plt.title("Mean Volume for L1")
    plt.xlabel("Index")
    plt.ylabel("Mean Volume")
    plt.legend()
    plt.grid(True)
    
    # プロットを保存
    plt.savefig(output_path)
    print(f"Plot saved to {output_path}")
    plt.close()

if __name__ == '__main__':
    # ファイルパス
    file_path = "/home/user/SA-EDS/dataset/L1_data_int_second.pkl"
    output_path = "/home/user/SA-EDS/results/plot_same_l1.png"

    try:
        # ファイルを読み込む
        with open(file_path, "rb") as f:
            data = pickle.load(f)

        # dicの初期化とデータの整理
        dic = {'277': [], '308': [], '328': [], '358': []}
        max_val = 0
        max_key = ""
        for key in data:
            if max_val < data[key]["mean_volume"]:
                max_val = data[key]["mean_volume"]
                max_key = key
            dic[key[0]].append(data[key]["mean_volume"])

        print(max_val, max_key)
        # dicをプロットして保存
        plot_and_save(dic, output_path)

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
